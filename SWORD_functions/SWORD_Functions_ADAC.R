
# SWORD ---------------------------------------------------------------
# Fits an SVM Rotation Tree random forest (bagging ensemble of TORS trees)
# with optional out-of-bag (OOB) error estimation. Supports parallel execution
# via furrr/future and optional tree-level timeout to skip slow trees.
#
# Args: (see inline documentation in the function signature)
#   Covariates   : data.frame of predictors
#   y            : numeric response vector
#   nmin         : min node size to attempt a split (default 5)
#   minleaf      : min leaf size (default 2)
#   cp           : complexity parameter (default 0.0)
#   n_perc       : number of quantile thresholds per split (default 1)
#   n_topCor     : number of top-correlated features per split (default 2)
#   threshold_COR: min |correlation| to pass feature filter (default 1)
#   m            : number of trees (default 10)
#   rf_var       : features sub-sampled per split (default ceil(p/3))
#   rand_ntopcor : randomly draw top-feature cardinality (default TRUE)
#   relation     : "Pearson" or "MI" for feature selection
#   Weight_Scheme: "scale" or "robust" observation weighting
#   type_of_svm  : "C-classification" or "nu-classification"
#   cost_C       : SVM cost C (default 1)
#   cost_nu      : SVM nu (default 0.5)
#   seed_BS      : bootstrap seed offset (default 25)
#   parallel     : use furrr parallel backend (default FALSE)
#   chunk        : process trees in chunks (default FALSE)
#   n_chunks     : number of chunks (NULL = auto)
#   OOB          : compute OOB predictions and metrics (default TRUE)
#   n_workers    : parallel workers (default 1)
#   verbose      : print tree index during sequential run (default FALSE)
#   timeout_tree : max seconds per tree; NULL = no limit (default 300)
#
# Returns:
#   List with elements: trees, oob_predictions, MSE, RMSE, MAE, Rsquared,
#   OOB_matrix, oob_errors_per_iter, time_tree, time_oob, n_skipped
#   (OOB metrics only present when OOB = TRUE and at least one tree succeeded)

SWORD <- function(
    Covariates,              # Matrix/data.frame of predictors
    y,                       # Response vector
    nmin = 5,                # Minimum number of samples required to attempt a split
    minleaf = 2,             # Minimum number of samples allowed in each terminal node
    cp = 0.0,                # Complexity parameter 
    n_perc = 1,              # number of quantiles to be used in each split
    n_topCor = 2,            # Number of top correlated features to retain
    threshold_COR = 1,       # Correlation threshold to filter features (use always 1)
    m = 10,                  # Number of trees in the forest
    rf_var = ceiling(ncol(Covariates)/1), # Number of variables sampled when searching for splits
    rand_ntopcor = TRUE,    # Whether to randomly pick the cardinality of top correlated features or just use n_topCor
    relation = "Pearson",    # "MI",  Dependency measure used (Pearson, Mutual information)
    Weight_Scheme="scale",     #"robust",  Weighting Scheme in the Weighted-SVM 
    type_of_svm="C-classification",# "nu-classification", Type of SVM 
    cost_C=1,
    cost_nu=0.5,
    seed_BS = 25,            # Seed for bootstrap sampling (reproducibility)
    parallel = FALSE,        # Enable parallel computation
    chunk = FALSE,           # Whether to process trees in chunks
    n_chunks = NULL,         # Number of chunks to split computation into (if enabled)
    OOB = TRUE,              # Whether to compute out-of-bag predictions
    n_workers = 1,               # Number of parallel workers (default 1 to avoid nested parallelism)
    verbose = FALSE,         # Verbose mode
    timeout_tree = 300       # Max seconds per tree (default 5 min). NULL = no timeout
) {
  
  Covariates$y <- y
  n <- nrow(Covariates)
  
  index_list <- lapply(1:m, function(i) {
    set.seed(i + seed_BS)
    boot_idx <- sample(n, n, replace = TRUE)
    oob_idx <- setdiff(1:n, boot_idx)
    list(boot = boot_idx, oob = oob_idx, seed = i + seed_BS)
  })
  
  build_tree <- function(index_set) {
    set.seed(index_set$seed)
    
    baggingData <- Covariates[index_set$boot, ]
    oobData <- Covariates[index_set$oob, ]
    
    y_bagging <- baggingData$y
    baggingData$y <- NULL
    oobData$y <- NULL
    
    # ---- train time (with timeout) ----
    tempo_albero_start <- proc.time()
    
    tree <- tryCatch({
      if (!is.null(timeout_tree)) setTimeLimit(elapsed = timeout_tree, transient = TRUE)
      res <- TORS(
        baggingData, y_bagging,
        nmin = nmin, minleaf = minleaf, cp = cp,
        n_perc = n_perc, n_topCor = n_topCor,
        threshold_COR = threshold_COR, rf_var = rf_var,
        rand_ntopcor = rand_ntopcor, relation = relation,
        Weight_Scheme = Weight_Scheme,
        type_of_svm = type_of_svm,
        cost_C = cost_C, cost_nu = cost_nu
      )
      if (!is.null(timeout_tree)) setTimeLimit(elapsed = Inf, transient = FALSE)
      res
    }, error = function(e) {
      setTimeLimit(elapsed = Inf, transient = FALSE)
      if (grepl("time|elapsed", e$message, ignore.case = TRUE)) {
        message("  [TIMEOUT] Tree seed=", index_set$seed, " exceeded ", timeout_tree, "s - skipped")
      } else {
        message("  [ERROR] Tree seed=", index_set$seed, ": ", e$message)
      }
      return(NULL)
    })
    
    tempo_albero <- as.numeric((proc.time() - tempo_albero_start)[["elapsed"]])
    
    # On timeout/error, return a NULL marker
    if (is.null(tree)) {
      return(list(
        tree = NULL,
        oob_preds = rep(NA_real_, n),
        time_tree = tempo_albero,
        time_oob = 0,
        skipped = TRUE
      ))
    }
    
    # ---- oob time + preds
    tempo_oob <- NA_real_
    oob_preds <- rep(NA_real_, n)
    
    if (OOB && length(index_set$oob) > 0) {
      tempo_oob_start <- proc.time()
      pred <- predict_TORS(oobData, tree)
      tempo_oob <- as.numeric((proc.time() - tempo_oob_start)[["elapsed"]])
      
      oob_preds[index_set$oob] <- pred
    } else {
      tempo_oob <- 0
    }
    
    list(
      tree = tree,
      oob_preds = oob_preds,
      time_tree = tempo_albero,
      time_oob = tempo_oob,
      skipped = FALSE
    )
  }
  
  if (parallel) {
    # Plan already configured externally - do not override
    # future::plan(future::multisession, workers = n_workers)  # ❌ COMMENTATO
    
    if (chunk) {
      if (is.null(n_chunks)) n_chunks <- ceiling(m / n_workers)
      index_chunks <- split(index_list, ceiling(seq_along(index_list) / ceiling(m / n_chunks)))
      
      results <- list()
      progressr::with_progress({
        for (ch in index_chunks) {
          tmp_results <- furrr::future_map(
            ch, 
            build_tree,
            .options = furrr::furrr_options(
              seed = TRUE,
              globals = TRUE,
              scheduling = Inf,
              chunk_size = NULL
            ),
            .progress = TRUE
          )
          results <- c(results, tmp_results)
          gc()
        }
      })
    } else {
      progressr::with_progress({
        results <- furrr::future_map(
          index_list, 
          build_tree,
          .options = furrr::furrr_options(
            seed = TRUE,
            globals = TRUE,
            scheduling = Inf,
            chunk_size = NULL
          ),
          .progress = TRUE
        )
      })
    }
    
    # Filter skipped trees (timeout/error)
    skipped <- vapply(results, function(r) isTRUE(r$skipped), logical(1))
    n_skipped <- sum(skipped)
    if (n_skipped > 0) message("  Skipped trees (timeout): ", n_skipped, " / ", m)
    valid <- results[!skipped]
    
    treeList <- lapply(valid, `[[`, "tree")
    
    # per-tree times (all, including skipped)
    time_tree <- vapply(results, `[[`, numeric(1), "time_tree")
    time_oob  <- vapply(results, `[[`, numeric(1), "time_oob")
    
    if (OOB && length(valid) > 0) {
      OOB_matrix <- sapply(valid, `[[`, "oob_preds")
      oob_counts <- rowSums(!is.na(OOB_matrix))
      oob_sums <- rowSums(OOB_matrix, na.rm = TRUE)
      final_oob_predictions <- oob_sums / oob_counts
      
      m_valid <- length(valid)
      oob_errors_per_iter <- sapply(1:m_valid, function(i) {
        partial <- OOB_matrix[, 1:i, drop = FALSE]
        row_mean <- rowSums(partial, na.rm = TRUE) / rowSums(!is.na(partial))
        mean((row_mean - y)^2, na.rm = TRUE)
      })
      
      mse_oob <- mean((final_oob_predictions - y)^2, na.rm = TRUE)
      rmse_oob <- sqrt(mse_oob)
      mae_oob <- mean(abs(final_oob_predictions - y), na.rm = TRUE)
      rsq_oob <- cor(final_oob_predictions, y, use = "complete.obs")^2
      
      return(list(
        trees = treeList,
        oob_predictions = final_oob_predictions,
        MSE = mse_oob,
        RMSE = rmse_oob,
        MAE = mae_oob,
        Rsquared = rsq_oob,
        OOB_matrix = OOB_matrix,
        oob_errors_per_iter = oob_errors_per_iter,
        time_tree = time_tree,
        time_oob = time_oob,
        n_skipped = n_skipped
      ))
    } else {
      return(list(
        trees = treeList,
        time_tree = time_tree,
        time_oob = time_oob,
        n_skipped = n_skipped
      ))
    }
    
  } else {
    treeList <- vector("list", m)
    
    # per-tree times (vectors)
    time_tree <- numeric(m)
    time_oob  <- numeric(m)
    n_skipped <- 0L
    
    oob_predictions <- rep(NA_real_, n)
    oob_counts <- rep(0L, n)
    oob_errors_per_iter <- rep(NA_real_, m)
    OOB_matrix <- matrix(NA_real_, nrow = n, ncol = m)
    
    for (i in 1:m) {
      if (verbose) message(i)
      
      set.seed(i + seed_BS)
      boot_idx <- sample(n, n, replace = TRUE)
      oob_idx <- setdiff(1:n, boot_idx)
      
      baggingData <- Covariates[boot_idx, ]
      oobData <- Covariates[oob_idx, ]
      
      y_bagging <- baggingData$y
      baggingData$y <- NULL
      oobData$y <- NULL
      
      # ---- train time (with timeout) ----
      tempo_albero_start <- proc.time()
      
      tree <- tryCatch({
        if (!is.null(timeout_tree)) setTimeLimit(elapsed = timeout_tree, transient = TRUE)
        res <- TORS(
          baggingData, y_bagging,
          nmin = nmin, minleaf = minleaf, cp = cp,
          n_perc = n_perc, n_topCor = n_topCor,
          threshold_COR = threshold_COR, rf_var = rf_var,
          rand_ntopcor = rand_ntopcor, relation = relation,
          Weight_Scheme = Weight_Scheme,
          type_of_svm = type_of_svm,
          cost_C = cost_C, cost_nu = cost_nu
        )
        if (!is.null(timeout_tree)) setTimeLimit(elapsed = Inf, transient = FALSE)
        res
      }, error = function(e) {
        setTimeLimit(elapsed = Inf, transient = FALSE)
        if (grepl("time|elapsed", e$message, ignore.case = TRUE)) {
          message("  [TIMEOUT] Tree ", i, " (seed=", i + seed_BS, ") exceeded ", timeout_tree, "s - skipped")
        } else {
          message("  [ERROR] Tree ", i, ": ", e$message)
        }
        return(NULL)
      })
      
      time_tree[i] <- as.numeric((proc.time() - tempo_albero_start)[["elapsed"]])
      
      # On timeout/error, skip this tree
      if (is.null(tree)) {
        n_skipped <- n_skipped + 1L
        next
      }
      
      treeList[[i]] <- tree
      
      # ---- oob time + preds
      if (OOB && length(oob_idx) > 0) {
        tempo_oob_start <- proc.time()
        oob_preds <- predict_TORS(oobData, tree)
        time_oob[i] <- as.numeric((proc.time() - tempo_oob_start)[["elapsed"]])
        
        oob_predictions[oob_idx] <- ifelse(
          is.na(oob_predictions[oob_idx]),
          oob_preds,
          oob_predictions[oob_idx] + oob_preds
        )
        oob_counts[oob_idx] <- oob_counts[oob_idx] + 1L
        OOB_matrix[oob_idx, i] <- oob_preds
        
        current_oob_preds <- oob_predictions / oob_counts
        oob_errors_per_iter[i] <- mean((current_oob_preds - y)^2, na.rm = TRUE)
      } else {
        time_oob[i] <- 0
      }
    }
    
    if (n_skipped > 0) message("  Skipped trees (timeout): ", n_skipped, " / ", m)
    
    # Remove NULL slots from treeList
    treeList <- Filter(Negate(is.null), treeList)
    
    # Remove all-NA columns from OOB_matrix (skipped trees)
    valid_cols <- colSums(!is.na(OOB_matrix)) > 0
    OOB_matrix <- OOB_matrix[, valid_cols, drop = FALSE]
    
    if (OOB && ncol(OOB_matrix) > 0) {
      final_oob_predictions <- oob_predictions / oob_counts
      
      m_valid <- ncol(OOB_matrix)
      oob_errors_per_iter <- sapply(1:m_valid, function(i) {
        partial <- OOB_matrix[, 1:i, drop = FALSE]
        row_mean <- rowSums(partial, na.rm = TRUE) / rowSums(!is.na(partial))
        mean((row_mean - y)^2, na.rm = TRUE)
      })
      
      mse_oob <- mean((final_oob_predictions - y)^2, na.rm = TRUE)
      rmse_oob <- sqrt(mse_oob)
      mae_oob <- mean(abs(final_oob_predictions - y), na.rm = TRUE)
      rsq_oob <- cor(final_oob_predictions, y, use = "complete.obs")^2
      
      return(list(
        trees = treeList,
        oob_predictions = final_oob_predictions,
        MSE = mse_oob,
        RMSE = rmse_oob,
        MAE = mae_oob,
        Rsquared = rsq_oob,
        OOB_matrix = OOB_matrix,
        oob_errors_per_iter = oob_errors_per_iter,
        time_tree = time_tree,
        time_oob = time_oob,
        n_skipped = n_skipped
      ))
    } else {
      return(list(
        trees = treeList,
        time_tree = time_tree,
        time_oob = time_oob,
        n_skipped = n_skipped
      ))
    }
  }
}




# get_deviance -----------------------------------------------------------------
# Computes the deviance (total sum of squares) of a numeric vector y.
# Equivalent to var(y) * (n - 1). Returns 0 for vectors of length <= 1.
#
# Args:
#   y : numeric vector
#
# Returns:
#   Non-negative numeric scalar
get_deviance <- function(y) {
  if (length(y) <= 1) return(0)
  
  # Calculate deviance using the formula: var(y) * (n - 1)
  deviance <- var(y) * (length(y) - 1)
  
  return(deviance)
}



# select_top_n_SVM -------------------------------------------------------------
# Returns the original indices of the n_top features with the highest scores,
# ranked in descending order.
#
# Args:
#   scores : numeric vector of feature scores
#   n_top  : integer, number of top features to select
#
# Returns:
#   Integer vector of length n_top with the original indices of top features
select_top_n_SVM <- function(scores, n_top) {
  
  # Create a data frame with scores and corresponding indices
  d <- data.frame(
    score  = data.table::copy(scores),
    index = seq(1, length(scores))
  )
  
  # Convert to data.table and order by descending scores
  data.table::setDT(d)
  data.table::setorder(d, -score)
  
  # Return the indices of the top n features
  top_n_indices <- d$index[1:n_top]
  
  return(top_n_indices)
}



# coef_unscaled_SVM ------------------------------------------------------------
# Extracts the SVM hyperplane coefficients and back-transforms them from the
# internally scaled feature space to the original (unscaled) feature space.
# Handles single-feature and multi-feature cases, and the case where the SVM
# was fitted without scaling.
#
# Args:
#   tree_SVM   : a fitted WeightSVM/e1071 SVM object (linear kernel)
#   Covariates : data.frame of training covariates used to fit tree_SVM
#
# Returns:
#   Named numeric vector of unscaled coefficients; the last element is "Int"
#   (the intercept)
coef_unscaled_SVM <- function(tree_SVM, Covariates) {
  
  # Get slope and intercept from scaled space
  w <- t(tree_SVM$coefs) %*% tree_SVM$SV
  intercept <- tree_SVM$rho
  
  if(dim(w)[2]==1){
    names(w)=paste(tree_SVM[["terms"]][[3]])
  
  # wname<<-w
  # Transform slope and intercept from scaled to unscaled space
  
  # Unscaled w and intercept
  m <- rep(0, ncol(Covariates))
  names(m) <- names(Covariates)  # Save variable names
  
  m[names(w)] <- w
  m_new <- as.matrix(m) / apply(Covariates, 2, sd,na.rm=T)
  # print((as.matrix(w)))
  # print(colMeans(Covariates[,colnames(w)]))
  # print(apply(Covariates[,colnames(w)], 2, sd,na.rm=T))
  
  int_new <- (-intercept - (as.matrix(w)) %*% 
                as.matrix(colMeans(Covariates[,names(w),drop=F]) / apply(Covariates[,names(w),drop=F], 2, sd,na.rm=T)))
  
  w_unscaled <- t(m_new)
  int_unscaled <- -as.numeric(int_new)
  
  coefficients <- c(w_unscaled, int_unscaled)
  coefficients[is.na(coefficients)]<-0
  names(coefficients) <- c(colnames(w_unscaled), "Int")
  }else{# Transform slope and intercept from scaled to unscaled space
    
    if(all(tree_SVM[["scaled"]] == FALSE)){
      m <- rep(0, ncol(Covariates))
      names(m) <- names(Covariates)  # Save variable names
      m[colnames(w)] <- w
      m_new=m
      int_new=intercept
      w_unscaled <- t(m_new)
      int_unscaled <- -as.numeric(int_new)
      
      coefficients <- c(w_unscaled, int_unscaled)
      coefficients[is.na(coefficients)]<-0
      names(coefficients) <- c(colnames(w_unscaled), "Int")
    }else{
    
    
    # Unscaled w and intercept
    m <- rep(0, ncol(Covariates))
    names(m) <- names(Covariates)  # Save variable names
    
    m[colnames(w)] <- w
    m_new <- as.matrix(m) / apply(Covariates, 2, sd,na.rm=T)
    
    int_new <- (-intercept - (as.matrix(w)) %*% 
                  as.matrix(colMeans(Covariates[,colnames(w)]) / apply(Covariates[,colnames(w)], 2, sd,na.rm=T)))
    
    w_unscaled <- t(m_new)
    int_unscaled <- -as.numeric(int_new)
    
    coefficients <- c(w_unscaled, int_unscaled)
    coefficients[is.na(coefficients)]<-0
    names(coefficients) <- c(colnames(w_unscaled), "Int")
    
    }
    
  }
  
  
  
  return(coefficients)
}



# best_split_SVM_ROT -----------------------------------------------------------
# Finds the best weighted linear SVM split for a node by:
#   1. Optionally sub-sampling features (random-forest style, via rf_var).
#   2. Filtering out constant features.
#   3. Selecting the most relevant features via Pearson correlation or Mutual
#      Information with the response y.
#   4. Scanning candidate dichotomization thresholds of y (n_perc quantiles)
#      and fitting a weighted SVM at each; with n_perc = 1 only the median is used.
#   5. Returning the SVM fitted at the threshold that minimises node deviance.
#
# Args:
#   x             : data.frame of covariates
#   y             : numeric response vector
#   n_perc        : number of candidate quantile thresholds (1 = median only)
#   n_topCor      : number of top-correlated features to retain
#   threshold_COR : minimum |correlation| for a feature to pass the filter
#   rf_var        : number of features to randomly sub-sample (NULL = all)
#   rand_ntopcor  : if TRUE, randomly draw the cardinality of top features
#   relation      : feature-selection criterion: "Pearson" or "MI"
#   Weight_Scheme : observation weights: "scale" (standardised) or "robust" (MAD)
#   type_of_svm   : "C-classification" or "nu-classification"
#   cost_C        : cost parameter C for C-SVM
#   cost_nu       : nu parameter for nu-SVM
#
# Returns:
#   A fitted WeightSVM object (linear kernel) at the optimal threshold
best_split_SVM_ROT <- function(x, y, n_perc, n_topCor,threshold_COR,
                               rf_var,rand_ntopcor,relation,
                               Weight_Scheme = c("scale","robust"),
                               type_of_svm   = c("C-classification","nu-classification"),
                               cost_C = 1,
                               cost_nu = 0.5) {
  
  if(!is.null(rf_var) && rf_var>=n_topCor){
    
    selected_columns <- sample(1:ncol(x), ceiling(rf_var)) # Indices of the selected columns
    x <- x[, selected_columns] # New dataset with selected columns
  }
  # Eliminate features that are constant (sd = 0)
  non_constant_features <- x[, which(apply(x, 2, sd) != 0),drop=FALSE]
  
  cors=data.frame(abs(t(cor(non_constant_features,y))))
  cors=cors[,cors >= threshold_COR,drop=F]
  
  cor_threshold_features=(non_constant_features[,names(cors),drop=FALSE])
  if (ncol(cor_threshold_features) < 1){
    
    if (is.null(dim(non_constant_features)) || ncol(non_constant_features) == 0) {
      # No variable selection needed if all features are constant or only one feature is non-constant
      selected_features <- x
    } else {
      
      if(relation=="Pearson"){
        if((rand_ntopcor==FALSE)){
          # Select the top features with the highest absolute correlation with y
          selected_features <- data.frame(
            non_constant_features[, 
                                  select_top_n_SVM(as.vector(
                                    abs(cor(y, non_constant_features))), 
                                    min(n_topCor, ncol(non_constant_features))),drop=FALSE])
        }else{
          n_topCor_rand=sample(seq(2,rf_var,1),size=1)
          
          selected_features <- data.frame(
            non_constant_features[, 
                                  select_top_n_SVM(as.vector(
                                    abs(cor(y, non_constant_features))), 
                                    min(n_topCor_rand, ncol(non_constant_features))),drop=FALSE])
          
        }
      }
      if(relation=="MI"){
        if((rand_ntopcor==FALSE)){
          # Select the top features with the highest absolute correlation with y
          discretized_features <- infotheo::discretize(non_constant_features)
          discretized_y <- infotheo::discretize(data.frame(y))  # y must be a vector or a column
          
          mi_scores <- sapply(seq_len(ncol(discretized_features)), function(i) {
            infotheo::mutinformation(discretized_features[, i, drop = FALSE], discretized_y)
          })
          
          # Feature selection with highest mutual information
          selected_features <- data.frame(
            non_constant_features[, 
                                  select_top_n_SVM(as.vector(mi_scores), 
                                                   min(n_topCor, ncol(non_constant_features))), 
                                  drop = FALSE]
          )
        }else{
          n_topCor_rand=sample(seq(2,rf_var,1),size=1)
          
          discretized_features <- infotheo::discretize(non_constant_features)
          discretized_y <- infotheo::discretize(data.frame(y))  # y must be a vector or a single column
          
          mi_scores <- sapply(seq_len(ncol(discretized_features)), function(i) {
            infotheo::mutinformation(discretized_features[, i, drop = FALSE], discretized_y)
          })
          
          # Feature selection with highest mutual information
          selected_features <- data.frame(
            non_constant_features[, 
                                  select_top_n_SVM(as.vector(mi_scores), 
                                                   min(n_topCor_rand, ncol(non_constant_features))), 
                                  drop = FALSE]
          )
          
        }
      }
      
    }
    
  }else{
    
    selected_features<-cor_threshold_features
  }
  
  if(n_perc==1){#no need to find best hyperplane, it is just one
    
    t=median(y,na.rm=TRUE)
    # Create the final SVM model with the optimal threshold
    data_svm_best <- data.frame(y_t = factor(y <= t), selected_features)
    if (length(unique(data_svm_best$y_t))==1){
      #pass from <= to just < so the dichotomization does not trasform it
      #in a one modality variable
      data_svm_best <- data.frame(y_t = factor(y < t), selected_features)
      
    }
    formula_SVM_Tree <- as.formula(paste("y_t ~ .", collapse = "+"))
    # Weight_Scheme: "scale" or "robust"
    
    if (Weight_Scheme == "scale") {
      weights <- abs(scale(y))
      
    } else if (Weight_Scheme == "robust") {
      
      med_y <- median(y)
      
      s <- mad(y)
      if (!is.finite(s) || s == 0) s <- sd(y)
      if (!is.finite(s) || s == 0) s <- 1
      
      weights <- abs(y - med_y) / s
      weights[weights == 0] <- 1e-10
      weights <- as.matrix(weights)
    } else {
      stop("Weight_Scheme must be 'scale' or 'robust'")
    }
    
    tree_SVM <- WeightSVM::wsvm(
      formula_SVM_Tree,
      data = data_svm_best,
      weight = weights,
      kernel = "linear",
      type = type_of_svm,
      cost=cost_C,
      nu = cost_nu,
      fitted = FALSE,
      scale = TRUE
    )
    
  }else{
    # Compute quantiles for splitting
    thr <- unique(quantile(y, seq(0, 1, max(1 / (n_perc + 1), (1 / (length(y) + 1)))), na.rm = TRUE))
    
    # Create an empty dataframe to store results
    out_svm <- data.frame(
      threshold = thr[-length(thr)],
      deviance = NA_real_
    )
    
    # SVM works with categorical variables, so split at different quantiles
    # and find the best split that minimizes deviance
    data_svm_list <- sapply(thr[-length(thr)], function(t) {
      
      # Dichotomize variable y
      data_svm <- data.frame(y_t = factor(y <= t), selected_features)
      
      # Create SVM formula
      svm_formula <- as.formula(paste("y_t ~ .", collapse = "+"))
      
      
      if (Weight_Scheme == "scale") {
        weights <- abs(scale(y))
        
      }else if (Weight_Scheme == "robust") {
        
        med_y <- median(y)
        
        s <- mad(y)
        if (!is.finite(s) || s == 0) s <- sd(y)
        if (!is.finite(s) || s == 0) s <- 1
        
        weights <- abs(y - med_y) / s
        weights[weights == 0] <- 1e-10
        weights <- as.matrix(weights)
      } else {
        stop("Weight_Scheme must be 'scale' or 'robust'")
      }
      
      
      # Train SVM model
      svmtree <- WeightSVM::wsvm(
        svm_formula,
        data = data_svm,
        weight = weights,
        kernel = "linear",
        type = type_of_svm,
        cost=cost_C,
        nu = cost_nu,
        fitted = FALSE,
        scale = TRUE
      )
      
      
      svm_split_side <- function(svmtree, X) {
        # X: data.frame of features used in the svm (selected_features)
        
        coefs <- coef_unscaled_SVM(svmtree, Covariates = X)
        
        b <- unname(coefs["Int"])
        w <- coefs[names(coefs) != "Int"]
        
        # align columns (in case some are missing)
        common <- intersect(names(w), colnames(X))
        if (length(common) == 0) {
          # degenerate split: assign all observations to the left branch
          return(rep(FALSE, nrow(X)))
        }
        
        w_use <- w[common]
        Xmat  <- as.matrix(X[, common, drop = FALSE])
        
        score <- drop(Xmat %*% w_use + b)
        
        # "TRUE" side = score >= 0 (convention)
        score >= 0
      }
      idx_right <- svm_split_side(svmtree, selected_features)
      
      deviances <- get_deviance(y[idx_right]) +
        get_deviance(y[!idx_right])
      
      
      return(deviances)
    })
    
    # Save deviances to the out_svm dataframe
    out_svm$deviance <- data_svm_list
    
    # Find the best deviance and corresponding threshold
    t_best <- out_svm[which(out_svm$deviance == min(out_svm$deviance, na.rm = TRUE))[1], "threshold"]
    dev_SVM_best <- min(out_svm$deviance, na.rm = TRUE)
    
    # Create the final SVM model with the optimal threshold
    data_svm_best <- data.frame(y_t = factor(y <= t_best), selected_features)
    formula_SVM_Tree <- as.formula(paste("y_t ~ .", collapse = "+"))
    
    if (Weight_Scheme == "scale") {
      weights <- abs(scale(y))
      
    } else if (Weight_Scheme == "robust") {
      
      med_y <- median(y)
      
      s <- mad(y)
      if (!is.finite(s) || s == 0) s <- sd(y)
      if (!is.finite(s) || s == 0) s <- 1
      
      weights <- abs(y - med_y) / s
      weights[weights == 0] <- 1e-10
      weights <- as.matrix(weights)
    } else {
      stop("Weight_Scheme must be 'scale' or 'robust'")
    }
    
    
    tree_SVM <- WeightSVM::wsvm(
      formula_SVM_Tree,
      data = data_svm_best,
      weight = weights,
      kernel = "linear",
      type = type_of_svm,
      cost = cost_C,
      nu = cost_nu,
      fitted = FALSE,
      scale = TRUE
    )
  }
  # The result of the function is the SVM classifier
  return(tree_SVM)
}

# node_split_SVM_ROT -----------------------------------------------------------
# Applies best_split_SVM_ROT to partition a node's data into left and right
# branches according to the sign of the SVM decision function.
#
# Args:
#   Covariates : data.frame of covariates
#   response   : numeric response vector
#   (remaining arguments are forwarded to best_split_SVM_ROT)
#
# Returns:
#   List with elements:
#     left_df  – data.frame for the left  branch (response in last column)
#     right_df – data.frame for the right branch (response in last column)
#     n_left   – number of observations in the left  branch
#     n_right  – number of observations in the right branch
#     tree_SVM – the fitted SVM object used for the split
node_split_SVM_ROT <- function(Covariates, response, n_perc = n_perc,
                               n_topCor = n_topCor,threshold_COR = threshold_COR,rf_var=rf_var,
                               rand_ntopcor = rand_ntopcor,relation=relation,
                               Weight_Scheme=Weight_Scheme,
                               type_of_svm=type_of_svm,cost_C=cost_C,cost_nu=cost_nu) {
  
  # Find the best SVM split
  tree_SVM <- best_split_SVM_ROT(x = Covariates, y = response, n_perc = n_perc,
                                 n_topCor = n_topCor, threshold_COR = threshold_COR,
                                 rf_var=rf_var,rand_ntopcor=rand_ntopcor,relation=relation,
                                 Weight_Scheme=Weight_Scheme,type_of_svm=type_of_svm,cost_C=cost_C,cost_nu=cost_nu)
  
  # # Merge Covariates and response
  merged_data <- data.frame(Covariates, y = response)
  # 

  
  predizioni=Pred_coef_SVM(Covariates,coef_unscaled_SVM(tree_SVM,Covariates))
  left_data <- merged_data[predizioni == TRUE, ]
  right_data <- merged_data[predizioni == FALSE, ]
  
  return(list(
    "left_df" = left_data,
    "right_df" = right_data,
    "n_left" = nrow(left_data),
    "n_right" = nrow(right_data),
    "tree_SVM" = tree_SVM
  ))
}




# tree_grow_SVM_ROT ------------------------------------------------------------
# Recursively grows a single SVM Rotation Tree (TORS). At each internal
# node the best hyperplane split is found via node_split_SVM_ROT; the function
# then recurses on left and right branches. Growth stops when any of the
# following conditions is met:
#   - All observations share the same y value.
#   - The node has fewer than nmin observations.
#   - Either child branch would have fewer than minleaf observations.
#   - The relative deviance reduction (cp_node) falls below cp.
#
# Args:
#   Covariates        : data.frame of covariates
#   y                 : numeric response vector
#   nmin              : minimum node size to attempt a split
#   minleaf           : minimum leaf size
#   cp                : complexity parameter (minimum relative deviance reduction)
#   n_perc            : number of candidate thresholds per split
#   n_topCor          : number of top features to use per split
#   threshold_COR     : minimum |correlation| threshold for feature selection
#   original_deviance : root deviance, used to compute relative cp values
#   depth             : current depth counter (incremented at each recursion)
#   (remaining arguments forwarded to node_split_SVM_ROT)
#
# Returns:
#   Nested list representing a binary tree:
#     Internal node: list(Info = <node metadata>, Left_branch, Right_branch)
#     Leaf node:     list(mean, n, dev, depth, Leaf)
tree_grow_SVM_ROT <- function(
    Covariates, y, nmin = 20, minleaf = round(nmin / 3),
    cp = 0.01, n_perc = 100, n_topCor = ncol(Covariates),
    threshold_COR = threshold_COR,
    original_deviance = original_deviance, depth = 0,rf_var=rf_var,
    rand_ntopcor=rand_ntopcor,
    relation=relation,Weight_Scheme=Weight_Scheme,type_of_svm=type_of_svm,
    cost_C=cost_C,cost_nu=cost_nu) {
  
  depth = depth
  
  # If the node is pure or smaller than nmin, create a leaf node
  if (length(unique(y)) == 1 || nrow(Covariates) < nmin) {
    leaf = list(mean(y), length(y), get_deviance(y), depth, Leaf)
    names(leaf) = c("mean", "n", "dev", "depth", "Leaf")
    Leaf <<- Leaf + 1
    return("leaf" = leaf)
  }
  
  # Find the best SVM split
  nodes = node_split_SVM_ROT (Covariates, y, n_perc = n_perc, n_topCor = n_topCor,
                              threshold_COR = threshold_COR,rf_var=rf_var,
                              rand_ntopcor=rand_ntopcor,relation=relation,
                              Weight_Scheme=Weight_Scheme,type_of_svm=type_of_svm,cost_C=cost_C,cost_nu=cost_nu)
  tree_SVM = nodes[["tree_SVM"]]

  
  # Calculate deviances for left and right branches
  deviances = get_deviance(nodes[["left_df"]][, ncol(nodes[["left_df"]])]) +
    get_deviance(nodes[["right_df"]][, ncol(nodes[["right_df"]])])
  
  dev_tot = get_deviance(y)

  
  CC_par = deviances / original_deviance
  
  # Information about the node
  coef_unscaled_ = coef_unscaled_SVM(tree_SVM, Covariates)
  scaled_coef=coef_scaled(tree_SVM, Covariates)
  infonode = list(mean(y), length(y), deviances,get_deviance(y), dev_tot, CC_par, depth, coef_unscaled_,scaled_coef)
  names(infonode) = c("mean", "n", "dev","dev_presplit", "pre_splitDev", "cp", "depth", "coeffs","Scaled_coeff")
  
  #print(CC_par)
  # If the quality is not sufficient, create a leaf node
  if (length(nodes[["left_df"]][, 1]) < minleaf ||
      length(nodes[["right_df"]][, 1]) < minleaf ||
      CC_par < cp) {
    leaf = list(mean(y), length(y), get_deviance(y), depth, Leaf)
    names(leaf) = c("mean", "n", "dev", "depth", "Leaf")
    Leaf <<- Leaf + 1
    return("leaf" = leaf)
  }
  
  # Recursion for the left branch
  leftBranch = tree_grow_SVM_ROT(
    Covariates = nodes[["left_df"]][-c(ncol(nodes[["left_df"]]) )],
    y = nodes[["left_df"]][, ncol(nodes[["left_df"]])],
    nmin, minleaf, cp, n_perc, n_topCor, threshold_COR
    , original_deviance, depth = depth + 1,rf_var=rf_var,
    rand_ntopcor=rand_ntopcor,relation=relation,Weight_Scheme=Weight_Scheme,
    type_of_svm =type_of_svm ,cost_C=cost_C,cost_nu=cost_nu
  )
  
  # Recursion for the right branch
  rightBranch = tree_grow_SVM_ROT(
    Covariates = nodes[["right_df"]][-c(ncol(nodes[["right_df"]]) )],
    y = nodes[["right_df"]][, ncol(nodes[["right_df"]])],
    nmin, minleaf, cp, n_perc, n_topCor, threshold_COR
    , original_deviance, depth = depth + 1,rf_var=rf_var,
    rand_ntopcor=rand_ntopcor,relation=relation,Weight_Scheme=Weight_Scheme,
    type_of_svm =type_of_svm,cost_C=cost_C,cost_nu=cost_nu
  )
  
  
  return(
    list( "Info" = infonode,
          "Left_branch" = c(leftBranch), "Right_branch" = c(rightBranch))
  )
}

# TORS ----------------------------------------------------------------------
# Fits a single SVM Rotation Tree for regression. Initialises the global Leaf
# counter, computes the root deviance, calls tree_grow_SVM_ROT, then cleans up.
#
# Args:
#   Covariates    : data.frame of covariates
#   y             : numeric response vector
#   nmin          : minimum node size to attempt a split (default 5)
#   minleaf       : minimum leaf size (default round(nmin/3))
#   cp            : complexity parameter (default 0.01)
#   n_perc        : number of quantile thresholds per split (default 1)
#   n_topCor      : number of top-correlated features to select (default 2)
#   threshold_COR : minimum |correlation| threshold (default 1 = all features)
#   rf_var        : number of features to randomly sub-sample per split
#   rand_ntopcor  : randomly draw the top-feature cardinality (default FALSE)
#   relation      : feature-selection criterion: "Pearson" or "MI"
#   Weight_Scheme : observation weighting: "scale" or "robust"
#   type_of_svm   : "C-classification" or "nu-classification"
#   cost_C        : SVM cost parameter C
#   cost_nu       : SVM nu parameter
#
# Returns:
#   Nested list (binary tree) as produced by tree_grow_SVM_ROT
TORS <- function(Covariates, y, nmin = 5, minleaf = round(nmin / 3),
                    cp = 0.01, n_perc = 1, n_topCor = 2,
                    threshold_COR = 1,rf_var=ncol(Covariates),rand_ntopcor=FALSE,
                    relation="Pearson",Weight_Scheme="scale",
                    type_of_svm="C-classification",cost_C=1,cost_nu=0.5) {
  

  
  
  
  
  # Initialize Leaf counter
  Leaf <<- 1
  
  # Calculate the original deviance
  original_deviance = get_deviance(y)
  
  # Grow the SVM Rotation Tree
  svm_rot_tree = tree_grow_SVM_ROT(
    Covariates = Covariates,
    y = y,
    nmin = nmin,
    minleaf = minleaf,
    cp = cp,
    n_perc = n_perc,
    n_topCor = n_topCor,
    threshold_COR = threshold_COR,
    original_deviance = original_deviance,
    depth = 0,
    rf_var=rf_var,
    rand_ntopcor=rand_ntopcor,
    relation=relation,
    Weight_Scheme=Weight_Scheme,type_of_svm=type_of_svm,cost_C=cost_C,cost_nu=cost_nu)
  
  # Remove Leaf variable from the global environment
  rm(Leaf, pos = ".GlobalEnv")
  
  # Return the SVM Rotation Tree
  return(svm_rot_tree)
}


# Pred_coef_SVM ----------------------------------------------------------------
# Classifies observations using a linear hyperplane defined by a named
# coefficient vector (including intercept "Int"). Returns TRUE for observations
# on the negative side of the hyperplane (decision value < 0), FALSE otherwise.
# This convention determines left/right branch assignment within TORS.
#
# Args:
#   data   : data.frame of covariates (column names must match names(coeffs))
#   coeffs : named numeric vector including an element named "Int" (intercept)
#
# Returns:
#   Logical vector of length nrow(data)
Pred_coef_SVM <- function(data, coeffs) {
  coeffs=data.frame(t(coeffs))
  coefficients=coeffs[,-which(names(coeffs)=="Int")]
  intercept=coeffs[,which(names(coeffs)=="Int"),drop=F]
  # Check that the number of coefficients matches the columns of the dataset
  if (length(coefficients) != ncol(data)) {
    stop("The number of coefficients must match the number of columns in the dataset.")
  }
  
  # Compute the dot product between the coefficients and the data, adding the intercept
  
  decision_values <- (as.matrix(data)) %*% as.numeric(coefficients) - as.numeric(intercept)
  
  # Create a vector that identifies the side of the hyperplane
  
  group <- ifelse(decision_values < 0, TRUE, FALSE)
  
  # Add the group to the dataset
  
  data$Group <- group
  
  # Return the dataset split into two groups
  return(data$Group)
}


# ==== Fast prediction helpers ================================================

# precompute_predictions_SVM ---------------------------------------------------
# Recursively traverses an TORS tree and pre-computes predictions for all
# rows of data in a single pass (avoids per-row tree traversal overhead).
#
# Args:
#   tree   : TORS tree (nested list with Left_branch/Right_branch, or leaf)
#   data   : data.frame of covariates to predict on
#   output : name of the leaf field to return as prediction (default "mean")
#
# Returns:
#   Numeric vector of length nrow(data) with predicted values
precompute_predictions_SVM <- function(tree, data,output = "mean") {
  if (names(tree[1]) == "mean") {
    return(rep(tree[[output]], nrow(data)))
  } else {
    left_preds <- precompute_predictions_SVM(tree[["Left_branch"]], data,output)
    right_preds <- precompute_predictions_SVM(tree[["Right_branch"]], data,output)
    return(ifelse(Pred_coef_SVM(data, tree[["Info"]][["coeffs"]] ) == TRUE, left_preds, right_preds))
  }
}

# go_tree_svm_rot_optimized ----------------------------------------------------
# Retrieves the pre-computed prediction for a single observation index.
#
# Args:
#   row_number  : integer row index
#   predictions : numeric vector produced by precompute_predictions_SVM
#
# Returns:
#   Scalar numeric prediction for that observation
go_tree_svm_rot_optimized <- function(row_number, predictions) {
  return(predictions[row_number])
}

# predict_TORS ---------------------------------------------------
# Main prediction interface for a single TORS tree. Calls
# precompute_predictions_SVM for efficient batch prediction.
#
# Args:
#   covariates : data.frame of covariates to predict on
#   tree       : fitted TORS tree object
#   output     : leaf field to use as prediction (default "mean")
#
# Returns:
#   Numeric vector of predictions, one per row of covariates
predict_TORS <- function(covariates, tree, output = "mean") {
  predictions <- precompute_predictions_SVM(tree, covariates,output)
  return(predictions)
}



# stripname --------------------------------------------------------------------
# Recursively removes all list elements whose name matches `name` from a nested
# list (e.g., strips "Info" metadata nodes from an TORS tree before plotting).
#
# Args:
#   x    : nested list (e.g., an TORS tree)
#   name : character string, name of elements to remove at every level
#
# Returns:
#   Nested list with all matching named elements removed
stripname <- function(x, name) {
  this_depth <- depth_strip(x)
  
  if (this_depth == 0) {
    return(x)
  } else if (length(name_index <- which(names(x) == name))) {
    x <- x[-name_index]
  }
  
  return(lapply(x, stripname, name))
}


# depth_strip ------------------------------------------------------------------
# Computes the maximum nesting depth of a (possibly nested) list structure.
# Source: http://stackoverflow.com/questions/13432863
#
# Args:
#   this       : any R object (list or atomic value)
#   this_depth : internal recursion counter; start at 0
#
# Returns:
#   Integer: maximum depth of nested lists (0 for atomic values)
depth_strip <- function(this, this_depth = 0) {
  if (!is.list(this)) {
    return(this_depth)
  } else {
    return(max(unlist(lapply(this, depth_strip, this_depth = this_depth + 1))))
  }
}



# Depth_tree -------------------------------------------------------------------
# Returns the maximum depth of an TORS tree by converting it to a
# data.tree object and reading the "depth" field.
#
# Args:
#   svm_tree : fitted TORS tree object
#
# Returns:
#   Integer: maximum node depth across the tree
Depth_tree <- function(svm_tree) {
  
  pruned_svm <- stripname(svm_tree, "SPLIT")
  
  dt_svm <- data.tree::FromListSimple(pruned_svm)
  
  df_dt_svm <- as.data.frame(dt_svm, row.names = NULL, optional = FALSE, 
                             "mean", "cp", 'height', 'level', "mean", "cp", "n", "dev", "Leaf", "depth")
  
  return(max(df_dt_svm$depth, na.rm = TRUE))
}

# N_of_Leafs -------------------------------------------------------------------
# Returns the number of leaf nodes in an TORS tree.
#
# Args:
#   svm_tree : fitted TORS tree object
#
# Returns:
#   Integer: number of leaves
N_of_Leafs <- function(svm_tree) {
  
  pruned_svm <- stripname(svm_tree, "SPLIT")
  
  dt_svm <- data.tree::FromListSimple(pruned_svm)
  
  df_dt_svm <- as.data.frame(dt_svm, row.names = NULL, optional = FALSE, 
                             "mean", "cp", 'height', 'level', "mean", "cp", "n", "dev", "Leaf", "depth")
  
  return(max(df_dt_svm$Leaf, na.rm = TRUE))
}



# Leaf_or_Depth_for_CV ---------------------------------------------------------
# Extracts a tree complexity statistic (number of leaves or maximum depth) from
# a nested list of cross-validation results and returns it as a matrix.
#
# Args:
#   lista_risultati : nested list of CV results, indexed as [fold][[seed]][[type]]
#   stat            : "Leaf" (number of leaves) or "Depth" (maximum depth)
#   type_tree       : 1 = rpart tree, 2 = TORS tree
#   times_metrica   : number of outer CV repetitions (rows of the output matrix)
#
# Returns:
#   Numeric matrix [times_metrica x 10] of the requested complexity statistic
Leaf_or_Depth_for_CV <- function(lista_risultati,stat="Leaf",type_tree=type_tree,times_metrica=10) {
  matrix_result = matrix(NA, times_metrica, 10)
  
  for (i in 1:length(lista_risultati)) {
    for (k in 1:length(lista_risultati[[i]])) {
      
      
      if(stat=="Leaf" & type_tree==1 ){
        matrix_result[i, k]=  sum(lista_risultati[[i]][[k]][[type_tree]][["model"]]$frame$ncompete == 0)
      } else if (stat=="Depth" & type_tree==1){
        matrix_result[i, k]=  max(rpart:::tree.depth( as.numeric(
          rownames(lista_risultati[[i]][[k]][[type_tree]][["model"]]$frame))))
      }
      
      
      if(stat=="Leaf" & type_tree==2 ){
        matrix_result[i, k]=  N_of_Leafs(lista_risultati[[i]][[k]][[type_tree]][["model"]])
      } else if (stat=="Depth"& type_tree==2){
        matrix_result[i, k]=  Depth_tree(lista_risultati[[i]][[k]][[type_tree]][["model"]])
      }
      
      
    }
  }
  
  return(matrix_result)
}




# depth_svm_rot ----------------------------------------------------------------
# Recursively computes the nesting depth (number of levels) of an TORS tree
# or any nested list structure.
#
# Args:
#   x : nested list (TORS tree or any R list)
#
# Returns:
#   Integer: maximum nesting depth (1 + max depth of children)
depth_svm_rot <- function(x) {
  if (is.list(x)) {
    return(1 + max(sapply(x, depth_svm_rot)))
  } else {
    return(0)
  }
}


# toTree -----------------------------------------------------------------------
# Converts leaf nodes of an TORS tree into a named nested list suitable for
# data.tree rendering. Each leaf gets a descriptive label summarising mean, dev,
# and n; internal nodes are passed through recursively unchanged.
#
# Args:
#   x : TORS tree or subtree (nested list)
#
# Returns:
#   Nested list formatted for data.tree::FromListSimple
toTree <- function(x) {
  d <- depth_svm_rot(x)
  if(d > 1) {
    lapply(x, toTree)
  } else {
    children2 = (lapply((x), function(nm) list(value=nm)))
    nomi_child=list(
      paste("mean",round(children2[["mean"]][["value"]],2),
            paste("dev",round(children2[["dev"]][["value"]],2)),
            paste("n",round(children2[["n"]][["value"]],2)),sep=", "),
      
      paste("dev",round(children2[["dev"]][["value"]],2)),
      paste("n",round(children2[["n"]][["value"]],2)),
      #paste("Leaf",round(children2[["Leaf"]][["value"]],2)),
      paste("depth",round(children2[["depth"]][["value"]],2)) 
    )
    
    children=list()
    children[[nomi_child[[1]]]]=list("value"=paste("mean",round(children2[["mean"]][["value"]],2),
                                                   paste("dev",round(children2[["dev"]][["value"]],2)),
                                                   paste("n",round(children2[["n"]][["value"]],2))
                                                   
    ))
    
    
    children2=children
    
    
  }
}


# plot_SVM.ROT -----------------------------------------------------------------
# Converts an TORS tree into a data.tree Node object for visualisation.
# Strips the specified metadata field before conversion; set strip = NULL to
# retain internal-node information.
#
# Args:
#   Tree_SVM : fitted TORS tree
#   strip    : name of the field to strip (default "Info"; NULL = keep all)
#
# Returns:
#   data.tree Node object
plot_SVM.ROT=function(Tree_SVM,strip="Info"){
  
  prun_svm=(stripname(Tree_SVM, "SPLIT"))
  prun_svm=(stripname(prun_svm, strip))
  dt <- data.tree::FromListSimple(toTree(prun_svm))
  return(dt)
  
}


# Table_SVM.ROT ----------------------------------------------------------------
# Converts an TORS tree into a flat data.frame where each row corresponds
# to one node. The raw "coeffs" string column is split into one numeric column
# per covariate plus an intercept column "Int".
#
# Args:
#   tree_SVM   : fitted TORS tree
#   Covariates : data.frame of covariates (column names label the coefficient columns)
#
# Returns:
#   data.frame with columns: mean, cp, height, level, n, dev, pre_splitDev,
#   Leaf, depth, <covariate names>, Int
Table_SVM.ROT= function(tree_SVM,Covariates){
  
  prun_svm=(stripname(tree_SVM, "SPLIT"))
  
  dt_svm=(data.tree::FromListSimple(prun_svm))
  
  df_dt_svm=as.data.frame(dt_svm, row.names = NULL, optional = T, 
                          "mean", "cp",'height','level',"cp","n","dev","pre_splitDev","Leaf","depth","coeffs")
  
  df_dt_svm[,c(colnames(Covariates),"Int")]= stringr::str_split_fixed(df_dt_svm$coeffs, ",", ncol(Covariates)+1)
  df_dt_svm=df_dt_svm[,-which(colnames(df_dt_svm)=="coeffs")]
  return(df_dt_svm)
}




# predict_SWORD_all ----------------------------------------------------------
# Generates a prediction matrix from an TORS random forest, where each
# column contains the predictions of one tree. Supports optional parallel
# execution via parallel/doParallel.
#
# Args:
#   covariates : data.frame of covariates to predict on
#   treeList   : object returned by SWORD (must contain element $trees)
#   parallel   : logical, enable parallel computation (default FALSE)
#   n_cores    : number of cores for parallel execution
#
# Returns:
#   Numeric matrix [nrow(covariates) x m] of per-tree predictions
predict_SWORD_all <- function(covariates, treeList, parallel = FALSE, n_cores = parallel::detectCores() - 2) {
  treeList <- treeList[["trees"]]
  
  if (parallel) {
    # --- PARALLEL VERSION ---
    
    
    library(foreach)
    library(doParallel)
    
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Export the necessary functions to the workers
    
    clusterExport(cl, varlist = c("predict_TORS",
                                  "precompute_predictions_SVM",
                                  "Pred_coef_SVM"),
                  envir = environment())
    
    # Parallel loop
    pred_mat <- foreach(tree = treeList, .combine = cbind) %dopar% {
      predict_TORS(covariates, tree)
    }
    
    stopCluster(cl)
  } else {
    # --- SEQUENTIAL VERSION ---
    pred_mat <- sapply(treeList, function(tree) {
      predict_TORS(covariates, tree)
    })
  }
  
  return(pred_mat)
}




# predict_SWORD --------------------------------------------------------------
# Generates ensemble predictions from an TORS random forest by averaging
# per-tree predictions (row means of the prediction matrix).
#
# Args:
#   covariates : data.frame of covariates to predict on
#   treeList   : object returned by SWORD (must contain element $trees)
#
# Returns:
#   Numeric vector of length nrow(covariates) with averaged ensemble predictions
predict_SWORD <- function(covariates, treeList){
  treeList=treeList[["trees"]]
  pred_mat <- sapply(treeList, function(tree) {
    predict_TORS(covariates, tree)
  })
  pred_mat=rowMeans(pred_mat)
  return(pred_mat)
}



# plot_rmse_trees --------------------------------------------------------------
# Iteratively computes RMSE as trees are added one by one to the ensemble and
# produces a line plot of RMSE vs. number of trees.
#
# Args:
#   trees       : numeric matrix [n x m], each column = predictions of one tree
#   test_target : numeric vector of true response values (length n)
#
# Returns:
#   Called for its side-effect (base R line plot); invisibly returns NULL
plot_rmse_trees <- function(trees, test_target) {
  # trees: matrix of tree predictions (each column is a tree)
  # test_target: vector of true values for rmse computation

  # Iterative rmse computation
  rmse_values <- c()
  for (i in 1:ncol(trees)) {
    # Mean of predictions up to the i-th column
    preds <- rowMeans(trees[, 1:i, drop = FALSE])

    # Compute rmse
    rmse_value <-Metrics:: rmse(test_target, preds)

    # Save rmse
    rmse_values <- c(rmse_values, rmse_value)
  }

  # Plot rmse values
  plot(1:length(rmse_values), rmse_values, type = "o", pch = 16,
       xlab = "Number of Trees", ylab = "rmse",
       main = "rmse as a function of the number of trees")
}




# coef_scaled ------------------------------------------------------------------
# Extracts the SVM hyperplane coefficients in the internally scaled feature space
# (i.e., without back-transforming to original units). Contrast with
# coef_unscaled_SVM which returns coefficients in the original feature space.
# Handles single-feature, unscaled, and multi-feature SVM cases.
#
# Args:
#   tree_SVM   : a fitted WeightSVM/e1071 SVM object (linear kernel)
#   Covariates : data.frame of training covariates used to fit tree_SVM
#
# Returns:
#   Named numeric vector of scaled coefficients; the last element is "Int"
#   (negated rho from the SVM)
coef_scaled <- function(tree_SVM, Covariates) {
  
  # Get slope and intercept from scaled space
  w <- t(tree_SVM$coefs) %*% tree_SVM$SV
  intercept <- tree_SVM$rho
  
  if (dim(w)[2] == 1) {
    names(w) <- paste(tree_SVM[["terms"]][[3]])
    
    # Allocate coefficient vector
    m <- rep(0, ncol(Covariates))
    names(m) <- names(Covariates)
    m[names(w)] <- w
    
    # No division by sd: keep coefficients as they are
    w_unscaled <- t(as.matrix(m))
    int_unscaled <- -as.numeric(intercept)
    
    coefficients <- c(w_unscaled, int_unscaled)
    coefficients[is.na(coefficients)] <- 0
    names(coefficients) <- c(colnames(w_unscaled), "Int")
    
  } else {
    
    if (all(tree_SVM[["scaled"]] == FALSE)) {
      m <- rep(0, ncol(Covariates))
      names(m) <- names(Covariates)
      m[colnames(w)] <- w
      w_unscaled <- t(m)
      int_unscaled <- -as.numeric(intercept)
      
    } else {
      m <- rep(0, ncol(Covariates))
      names(m) <- names(Covariates)
      m[colnames(w)] <- w
      
      # Keep w and intercept unscaled
      w_unscaled <- t(as.matrix(m))
      int_unscaled <- -as.numeric(intercept)
    }
    
    coefficients <- c(w_unscaled, int_unscaled)
    coefficients[is.na(coefficients)] <- 0
    names(coefficients) <- c(colnames(w_unscaled), "Int")
  }
  
  return(coefficients)
}











# Table_SVM.ROT_scaled ---------------------------------------------------------
# Converts an TORS tree into a flat data.frame using scaled (internal-space)
# SVM coefficients. The "Scaled_coeff" comma-separated string column is parsed
# into one numeric column per covariate plus an intercept column "Int".
#
# Args:
#   tree_SVM : fitted TORS tree
#
# Returns:
#   data.frame with one row per node and one column per covariate plus "Int"
Table_SVM.ROT_scaled <- function(tree_SVM) {
  names_covariates=names(tree_SVM[["Info"]][["coeffs"]])[names(tree_SVM[["Info"]][["coeffs"]]) != "Int"]
  prun_svm <- stripname(tree_SVM, "SPLIT")
  dt_svm <- data.tree::FromListSimple(prun_svm)

  df <- as.data.frame(dt_svm, row.names = NULL,
                      "mean", "cp", "height", "level", "cp", "n",
                      "dev", "pre_splitDev", "Leaf", "depth", "Scaled_coeff")

  # Split Scaled_coeff only if needed (and only once)
  coeff_matrix <- do.call(rbind, strsplit(df$Scaled_coeff, ",", fixed = TRUE))
  colnames(coeff_matrix) <- c(names_covariates, "Int")

  df <- cbind(df[, !(colnames(df) %in% "Scaled_coeff")], coeff_matrix)
  return(df)
}

# variable_importance_SVM_ROT --------------------------------------------------
# Computes variable importance for a single TORS tree. For each internal
# node, the deviance reduction (pre_splitDev - dev) is distributed among
# features proportionally to their normalised absolute scaled coefficients.
# Importance scores are then summed across all nodes.
#
# Args:
#   tree_SVM : fitted TORS tree
#
# Returns:
#   Named numeric vector of variable importance scores, one per covariate
variable_importance_SVM_ROT <- function(tree_SVM) {
  names_covariates=names(tree_SVM[["Info"]][["coeffs"]])[names(tree_SVM[["Info"]][["coeffs"]]) != "Int"]
  tavola_scaled <- Table_SVM.ROT_scaled(tree_SVM)
  tavola_na <- tavola_scaled[!is.na(tavola_scaled$cp), ]

  coef_mat <- tavola_na[, names_covariates, drop = FALSE]
  coef_mat[] <- lapply(coef_mat, as.numeric)  # faster and more consistent
  abs_coef_mat <- abs(as.matrix(coef_mat))    # use matrix directly for speed

  row_sums <- rowSums(abs_coef_mat)
  row_sums[row_sums == 0] <- NA

  norm_coef_mat <- abs_coef_mat / row_sums
  norm_coef_mat[is.na(norm_coef_mat)] <- 0

  delta_dev <- as.numeric(tavola_na$pre_splitDev) - as.numeric(tavola_na$dev)
  gain_mat <- norm_coef_mat * delta_dev

  variable_importance <- colSums(gain_mat, na.rm = TRUE)
  names(variable_importance) <- names_covariates

  return(variable_importance)
}

# VI_SWORD ---------------------------------------------
# Computes forest-level variable importance by averaging per-tree importance
# scores from variable_importance_SVM_ROT across all valid trees. Supports
# parallel computation via future/future.apply.
#
# Args:
#   foresta  : object returned by SWORD (must contain element $trees)
#   parallel : logical, enable parallel computation (default FALSE)
#   workers  : number of parallel workers (default 1)
#
# Returns:
#   Named numeric vector of mean variable importance (one per covariate), or
#   NULL if no valid trees are available
VI_SWORD <- function(foresta, parallel = FALSE, workers = 1) {
  if (parallel) {
    # Set up the parallelization plan
    future::plan(future::multisession, workers = workers)
    var_imp_list <- future.apply::future_lapply(foresta$trees, function(tree) {
      tryCatch(
        variable_importance_SVM_ROT(tree),
        error = function(e) NULL
      )
    })
  } else {
    var_imp_list <- lapply(foresta$trees, function(tree) {
      tryCatch(
        variable_importance_SVM_ROT(tree),
        error = function(e) NULL
      )
    })
  }

  var_imp_list <- Filter(Negate(is.null), var_imp_list)
  if (length(var_imp_list) == 0) return(NULL)

  var_imp_df <- do.call(rbind, var_imp_list)
  colMeans(var_imp_df, na.rm = TRUE)
}


# ==== Helpers for plotting single tree ===============================================================

# Create the hyperplane string:  + a1·x1 - a2·x2 ... + Int = 0

# eq_string --------------------------------------------------------------------
# Formats a named coefficient vector as a human-readable hyperplane equation
# string of the form "a1·x1 - a2·x2 + ... + b = 0". Zero coefficients are
# dropped; optionally only the top_k largest terms are shown.
#
# Args:
#   betas             : named numeric vector (may include element "Int")
#   digits            : number of decimal places for rounding (default 3)
#   top_k             : if non-NULL, show only the top_k largest-magnitude terms
#   include_intercept : if TRUE, append the intercept term to the string
#
# Returns:
#   Character string representing the hyperplane equation
eq_string <- function(betas, digits = 3, top_k = NULL, include_intercept = TRUE) {
  # betas is a named numeric (may contain "Int")
  if (is.null(betas) || length(betas) == 0) return("No split (all coefficients = 0)")
  
  # separate the intercept if present
  b0 <- 0
  if ("Int" %in% names(betas)) {
    b0 <- betas[["Int"]]
    betas <- betas[names(betas) != "Int"]
  }
  
  # remove zero coefficients
  betas <- betas[abs(betas) > 0]
  if (!length(betas)) {
    return(if (include_intercept && b0 != 0)
      paste0(format(round(b0, digits), trim=TRUE), " = 0")
      else "No split (all coefficients = 0)")
  }
  
  # optional: keep only top_k for compactness
  if (!is.null(top_k) && top_k < length(betas)) {
    ord <- order(abs(betas), decreasing = TRUE)
    betas <- betas[ord][seq_len(top_k)]
  }
  
  terms <- paste0(
    ifelse(sign(betas) >= 0, "+ ", "- "),
    format(round(abs(betas), digits), trim = TRUE),
    "·", names(betas)
  )
  lhs <- sub("^\\+\\s", "", paste(terms, collapse = " "))
  
  if (include_intercept) {
    inter <- if (b0 == 0) "" else {
      paste0(" ", ifelse(b0 >= 0, "+ ", "- "), format(round(abs(b0), digits), trim = TRUE))
    }
    paste0(lhs, inter, " = 0")
  } else {
    paste0(lhs, " = c")
  }
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ==== Build list for data.tree ======================================

# build_tree_list --------------------------------------------------------------
# Recursively converts an TORS tree into a nested named list suitable for
# data.tree::FromListSimple. Internal nodes are labelled with sample size,
# impurity gain (Δimp), and the hyperplane equation; leaf nodes are labelled
# with mean, sample size, and deviance.
#
# Args:
#   node       : TORS tree node (internal or leaf)
#   use_scaled : if TRUE, use Scaled_coeff for the equation string (default TRUE)
#   digits     : decimal places in the equation string (default 3)
#   top_k      : maximum number of equation terms to display (NULL = all)
#
# Returns:
#   Named nested list ready for data.tree::FromListSimple
build_tree_list <- function(node, use_scaled = TRUE, digits = 3, top_k = NULL) {
  # Internal node if $Info exists, leaf otherwise
  if (!is.null(node$Info)) {
    # Get coefficients (scaled or raw)
    betas <- if (use_scaled && !is.null(node$Info$Scaled_coeff)) {
      node$Info$Scaled_coeff
    } else {
      node$Info$coeffs
    }
    
    # Build the hyperplane equation string
    eq <- eq_string(betas, digits = digits, top_k = top_k, include_intercept = TRUE)
    
    # Delta impurity (dev presplit - dev post)
    dev_pre  <- node$Info$dev_presplit %||% node$Info$pre_splitDev
    dev_post <- node$Info$dev
    d_imp <- if (!is.null(dev_pre) && !is.null(dev_post)) dev_pre - dev_post else NA
    
    n <- node$Info$n %||% NA
    
    label <- sprintf("Node | n=%s | Δimp=%s\n%s",
                     n,
                     ifelse(is.na(d_imp), "NA", format(round(d_imp, 3), big.mark=",")),
                     eq)
    
    # Children
    children <- list()
    if (!is.null(node$Left_branch))  children <- c(children, build_tree_list(node$Left_branch,  use_scaled, digits, top_k))
    if (!is.null(node$Right_branch)) children <- c(children, build_tree_list(node$Right_branch, use_scaled, digits, top_k))
    
    out <- list()
    out[[label]] <- children
    return(out)
    
  } else {
    # Leaf: label with mean, n, dev
    mean_val <- node$mean %||% NA
    n_val    <- node$n    %||% NA
    dev_val  <- node$dev  %||% NA
    lbl <- sprintf("Leaf | mean=%.3f | n=%s | dev=%s",
                   mean_val,
                   n_val,
                   ifelse(is.na(dev_val), "NA", format(round(dev_val, 3), big.mark=",")))
    # Leaf -> empty list as children
    out <- list()
    out[[lbl]] <- list()
    return(out)
  }
}

# ==== Wrapper plot ==========================================================

# plot_SVM_ROT -----------------------------------------------------------------
# Renders an interactive DiagrammeR visualisation of an TORS tree. Each
# internal node shows sample size, impurity gain, and the hyperplane equation;
# each leaf shows mean, sample size, and deviance. Uses a top-to-bottom layout.
#
# Args:
#   Tree_SVM   : fitted TORS tree
#   strip      : field to strip before plotting (default "Info")
#   use_scaled : if TRUE, display scaled-space coefficients (default TRUE)
#   digits     : decimal places in equation strings (default 3)
#   top_k      : maximum number of equation terms per node (default 6)
#
# Returns:
#   DiagrammeR graph object (rendered as a side-effect)
plot_SVM_ROT <- function(Tree_SVM, strip = "Info", use_scaled = TRUE, digits = 3, top_k = 6) {
  prun_svm <- Tree_SVM  # (optional stripname if needed)
  
  lst <- build_tree_list(prun_svm, use_scaled = use_scaled, digits = digits, top_k = top_k)
  
  library(data.tree)
  library(DiagrammeR)
  
  dt <- FromListSimple(lst)
  
  # Use direction = "climb" or "descend" (not "vertical")
  g <- ToDiagrammeRGraph(dt, direction = "climb")
  
  # Vertical orientation (Top→Bottom). For horizontal use "LR".
  g <- add_global_graph_attrs(
    graph = g,
    attr = "rankdir",
    value = "TB",   # "TB" (top-bottom), "LR" (left-right)
    attr_type = "graph"
  )
  
  DiagrammeR::render_graph(g)
}
