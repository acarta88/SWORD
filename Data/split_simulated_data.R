# ==============================================================================
# Split simulated_datasets_ADAC.rda into two parts (< 100 MB each for GitHub)
# Run this script ONCE from the SWORD/ folder before pushing to GitHub.
# ==============================================================================

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # set to Data/ folder

load("simulated_datasets_ADAC.rda")   # loads: list_simulated_df

seeds <- names(list_simulated_df)
cat("Total seeds:", length(seeds), "\n")
cat("Seeds:", paste(seeds, collapse = ", "), "\n\n")

mid <- ceiling(length(seeds) / 2)

list_simulated_df_part1 <- list_simulated_df[seeds[1:mid]]
list_simulated_df_part2 <- list_simulated_df[seeds[(mid + 1):length(seeds)]]

cat("Part 1 seeds:", paste(names(list_simulated_df_part1), collapse = ", "), "\n")
cat("Part 2 seeds:", paste(names(list_simulated_df_part2), collapse = ", "), "\n\n")

save(list_simulated_df_part1, file = "simulated_datasets_ADAC_part1.rda")
save(list_simulated_df_part2, file = "simulated_datasets_ADAC_part2.rda")

cat("Saved:\n")
cat("  simulated_datasets_ADAC_part1.rda :",
    round(file.size("simulated_datasets_ADAC_part1.rda") / 1e6, 1), "MB\n")
cat("  simulated_datasets_ADAC_part2.rda :",
    round(file.size("simulated_datasets_ADAC_part2.rda") / 1e6, 1), "MB\n")
