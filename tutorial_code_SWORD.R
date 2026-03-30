# ==============================================================================
# SWORD & TORS -- Tutorial Code
# Carta A., Frigau L.
# Paper: "Oblique Random Forests for Regression via Weighted Support Vector Machine"
# ADAC (Manuscript ID: ADAC-D-25-00591)
# ------------------------------------------------------------------------------
# SETUP: set the working directory to the root of the cloned repository, e.g.:
#   setwd("C:/path/to/SWORD")
# All paths below (Data/, SWORD_functions/) are relative to that root.
# ==============================================================================


# --- 1. Load packages and functions -------------------------------------------

library(data.table)# data.table      1.18.2.1
library(infotheo) # infotheo        1.2.0.1
library(WeightSVM)# WeightSVM       1.7.16
library(data.tree) # data.tree       1.2.0
library(randomForest) # randomForest    4.7.1.2
library(ggplot2)# ggplot2         4.0.2


source("SWORD_functions/SWORD_Functions_ADAC.R")

load("Data/simulated_datasets_ADAC.rda")

# Scenario: 500 obs, 10 features, 40% noise variables, 50% non-linearity
df <- list_simulated_df[["seed_1"]][["500_n_obs"]][["10_n_feat"]][["0.4_Noise"]][["0.5_NoLin"]][["0_cat"]][["1_err"]]

cat("Dataset:", nrow(df), "obs,", ncol(df) - 1, "features\n")


# --- 2. Train / test split (70 / 30) ------------------------------------------

y <- df[, 1];   X <- df[, -1, drop = FALSE]

set.seed(42)
train_idx <- sample(nrow(df), floor(0.7 * nrow(df)))
test_idx  <- setdiff(seq_len(nrow(df)), train_idx)

X_train <- X[train_idx, ]; y_train <- y[train_idx]
X_test  <- X[test_idx,  ]; y_test  <- y[test_idx]

cat("Train:", length(train_idx), "| Test:", length(test_idx), "\n")


# --- 3. Fit SWORD (forest, 100 trees) -----------------------------------------

set.seed(1)
sword_fit <- SWORD(
  Covariates    = X_train,
  y             = y_train,
  nmin          = 5,
  cp            = 0.0,
  n_perc        = 1,
  n_topCor      = 2,
  threshold_COR = 1,
  m             = 100,
  rf_var        = ncol(X_train),
  rand_ntopcor  = TRUE,
  relation      = "Pearson",
  Weight_Scheme = "scale",
  type_of_svm   = "C-classification",
  cost_C        = 1,
  OOB           = TRUE,
  parallel      = FALSE,
  verbose       = FALSE
)

pred_sword <- predict_SWORD(X_test, sword_fit)
rmse_sword <- sqrt(mean((pred_sword - y_test)^2))


# --- 4. Fit TORS (single tree) ------------------------------------------------
# TORS is deterministic when rand_ntopcor = FALSE, no additional seed needed.

tors_fit <- TORS(
  Covariates    = X_train,
  y             = y_train,
  nmin          = 5,
  cp            = 0.0,
  n_perc        = 1,
  n_topCor      = 2,
  threshold_COR = 1,
  rf_var        = ncol(X_train),
  rand_ntopcor  = FALSE,
  relation      = "Pearson",
  Weight_Scheme = "scale",
  type_of_svm   = "C-classification",
  cost_C        = 1
)

pred_tors <- predict_TORS(X_test, tors_fit)
rmse_tors <- sqrt(mean((pred_tors - y_test)^2))


# --- 5. Fit RF (benchmark) ----------------------------------------------------

set.seed(2)
rf_fit <- randomForest(
  x     = X_train,
  y     = y_train,
  ntree = 100,
  mtry  = ceiling(ncol(X_train) / 3)
)

pred_rf <- predict(rf_fit, X_test)
rmse_rf <- sqrt(mean((pred_rf - y_test)^2))


# --- 6. RMSE summary ----------------------------------------------------------

results <- data.frame(
  Model = c("SWORD", "TORS", "RF"),
  RMSE  = round(c(rmse_sword, rmse_tors, rmse_rf), 4)
)
print(results)


# --- 7. Variable Importance: SWORD vs RF --------------------------------------

vi_sword      <- VI_SWORD(sword_fit, parallel = FALSE)
vi_sword_norm <- vi_sword / sum(vi_sword)

vi_rf_norm    <- importance(rf_fit)[, "IncNodePurity"]
vi_rf_norm    <- vi_rf_norm / sum(vi_rf_norm)

make_vi_df <- function(vi_norm, model_label) {
  data.frame(
    variable = names(vi_norm),
    VI       = as.numeric(vi_norm),
    type     = ifelse(grepl("noise", names(vi_norm)), "Noise", "Informative"),
    model    = model_label
  )
}

vi_df          <- rbind(make_vi_df(vi_sword_norm, "SWORD"),
                        make_vi_df(vi_rf_norm,    "RF"))
var_order      <- names(sort(vi_sword_norm, decreasing = FALSE))
vi_df$variable <- factor(vi_df$variable, levels = var_order)

print(
  ggplot(vi_df, aes(x = variable, y = VI, fill = type, alpha = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Informative" = "#2166ac", "Noise" = "#d73027")) +
    scale_alpha_manual(values = c("SWORD" = 1, "RF" = 0.55)) +
    coord_flip() +
    facet_wrap(~ model, ncol = 2) +
    labs(
      title    = "Variable Importance: SWORD vs RF",
      subtitle = sprintf("RMSE -- SWORD: %.4f  |  TORS: %.4f  |  RF: %.4f",
                         rmse_sword, rmse_tors, rmse_rf),
      x = NULL, y = "Normalized VI",
      fill = "Variable type"
    ) +
    guides(alpha = "none") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
)


# --- 8. OOB diagnostics (SWORD) -----------------------------------------------

cat("OOB RMSE:", sword_fit$RMSE, "\n")

plot(sword_fit$oob_errors_per_iter, type = "l",
     xlab = "Number of trees", ylab = "OOB RMSE",
     main = "SWORD -- OOB convergence")
