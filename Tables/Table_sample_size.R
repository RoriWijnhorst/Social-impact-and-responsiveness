

# Define simulated ("true") values
sim_values <- data.frame(
  Parameter = c("B_0", "psi", "Sigma2_intercept", "Sigma2_psi", "Sigma2_phi", "cov_int_phi", "cov_int_psi", "cov_psi_phi"),
  Value = c(1, 0.3, 0.2, 0.1, 0.1, 0.08049845, -0.08485281, -0.07589466)
)

# List of result files and their labels
files <- c("Output/res_400_50_8_1x_M1.csv",
           "Output/res_800_100_8_1x_M1.csv",
           "Output/res_1600_200_8_1x_M1.csv",
           "Output/res_3200_400_8_1x_M1.csv",
           "Output/res_6400_800_8_1x_M1.csv")

labels <- c("400_50", "800_100", "1600_200", "3200_400",
            "6400_800")

# Initialize results table with true values
results_table <- sim_values

for (i in seq_along(files)) {
  df <- read.csv(files[i])
  
  # Compute mean posterior median per parameter
  means <- aggregate(X50. ~ Parameter, data = df, FUN = mean)
  
  # Merge with results table to get the true value for each parameter
  merged <- merge(results_table[, c("Parameter", "Value")], means, by = "Parameter", all.x = TRUE)
  
  # Calculate percentage bias
  merged$bias_percent <- 100 * (merged$X50. - merged$Value) / merged$Value
  
  # Add to the results table
  results_table <- merge(results_table, merged[, c("Parameter", "bias_percent")], by = "Parameter", all.x = TRUE)
  
  # Rename the new column
  colnames(results_table)[ncol(results_table)] <- labels[i]
}

# Optional: reorder rows to match the desired order from the example
param_order <-c("B_0", "psi", "Sigma2_intercept", "Sigma2_psi", "Sigma2_phi", "cov_int_phi", "cov_int_psi", "cov_psi_phi")


results_table <- results_table[match(param_order, results_table$Parameter), ]

# View the final table
print(results_table)
