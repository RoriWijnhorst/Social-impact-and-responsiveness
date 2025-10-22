library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(forcats)


df1_3200_1600_2_1x  <-read.csv("Output/res_3200_1600_2_1x_M1.csv") 
df1_3200_800_4_1x   <-read.csv("Output/res_3200_800_4_1x_M1.csv") 
df1_3200_400_8_1x   <-read.csv("Output/res_3200_400_8_1x_M1.csv") 
df1_3200_200_16_1x  <-read.csv("Output/res_3200_200_16_1x_M1.csv")
df1_3200_100_32_1x  <-read.csv("Output/res_3200_100_32_1x_M1.csv")

df2_3200_800_4_1x  <-read.csv("Output/res_3200_800_4_1x_M1.csv")
df2_3200_400_4_2x  <-read.csv("Output/res_3200_400_4_2x_M1.csv")
df2_3200_200_4_4x  <-read.csv("Output/res_3200_200_4_4x_M1.csv")
df2_3200_100_4_8x  <-read.csv("Output/res_3200_100_4_8x_M1.csv")

df3_3200_200_16_1x <-read.csv("Output/res_3200_200_16_1x_M1.csv")
df3_3200_200_8_2x  <-read.csv("Output/res_3200_200_8_2x_M1.csv")
df3_3200_200_4_4x  <-read.csv("Output/res_3200_200_4_4x_M1.csv")
df3_3200_200_2_8x  <-read.csv("Output/res_3200_200_2_8x_M1.csv")

df <- gdata::combine(
  df1_3200_1600_2_1x,  
  df1_3200_800_4_1x,  
  df1_3200_400_8_1x,  
  df1_3200_200_16_1x,
  df1_3200_100_32_1x,
  
  df2_3200_800_4_1x,
  df2_3200_400_4_2x,
  df2_3200_200_4_4x,
  df2_3200_100_4_8x,
  
  df3_3200_200_16_1x,
  df3_3200_200_8_2x,
  df3_3200_200_4_4x,
  df3_3200_200_2_8x)

# Means
Malpha   = 1      # mean focal effect
Mepsilon = 0      # mean partner effect
Mpsi     = 0.3    # interaction coefficient psi
Mx       = 0      # mean-centered covariate/opponent trait

# Variances
Valpha   = 0.2    # variance direct effect (mean behaviour)
Vepsilon = 0.01   # variance social partner effect (residual impact)
Vpsi     = 0.1    # variance in slopes (social responsiveness)
Vx       = 1      # variance social partner phenotype (impact covariate)

# Correlations between variance components
r_alpha_epsilon  =  0       # part of Cov int-phi
r_alpha_psi      =  -0.6    # part of Cov int-psi 
r_epsilon_psi    =  -0.6    # part of Cov psi-phi 
r_alpha_x        =  0.6     # part of Cov int-phi 
r_psi_x          =  -0.6    # part of Cov psi-phi 
r_epsilon_x      =  0

df$Q <- df$source
df$Q <- gsub('df','',df$Q)
df$Q <- as.factor(df$Q)
df[c('Q','N','ind','n_per_ind','repeats', 'Model')] <- str_split_fixed(df$Q, '_', 6)
df <- df %>%  mutate(ind = fct_relevel(ind,"1600", "800","400","200","100"))

Vphi <- Mpsi^2 * Vx + Vepsilon
Cov_alpha.phi <- r_alpha_epsilon*(sqrt(Valpha*Vepsilon))+ Mpsi*r_alpha_x*(sqrt(Valpha*Vx))
Cov_psi.phi <- r_epsilon_psi*(sqrt(Vepsilon*Vpsi)) + Mpsi*r_psi_x*(sqrt(Vpsi*Vx))
Cov_alpha.psi <- r_alpha_psi*(sqrt(Valpha*Vpsi))


df$RBias <- NA
df$RBias <- ifelse(df$Parameter=="psi", (df$X50.-Mpsi)/Mpsi, df$RBias)
df$RBias <- ifelse(df$Parameter=="Sigma2_x", (df$X50.-Vx)/Vx, df$RBias)
df$RBias <- ifelse(df$Parameter=="Sigma2_epsilon", (df$X50.-Vepsilon)/Vepsilon, df$RBias)
df$RBias <- ifelse(df$Parameter=="Sigma2_intercept", (df$X50.-Valpha)/Valpha, df$RBias)
df$RBias <- ifelse(df$Parameter=="Sigma2_psi", ((df$X50.-Vpsi)/Vpsi),df$RBias)
df$RBias <- ifelse(df$Parameter=="Sigma2_phi", ((df$X50.-Vphi)/Vphi),df$RBias)
df$RBias <- ifelse(df$Parameter=="cor1", r_alpha_epsilon,df$RBias)
df$RBias <- ifelse(df$Parameter=="cor2", ((df$X50.-r_alpha_psi)/r_alpha_psi),df$RBias)
df$RBias <- ifelse(df$Parameter=="cor3", ((df$X50.-r_epsilon_psi)/r_epsilon_psi),df$RBias)
df$RBias <- ifelse(df$Parameter=="cor4", ((df$X50.-r_alpha_x)/r_alpha_x),df$RBias)
df$RBias <- ifelse(df$Parameter=="cor5", ((df$X50.-r_psi_x)/r_psi_x),df$RBias)
df$RBias <- ifelse(df$Parameter=="cor6", r_epsilon_x,df$RBias)
df$RBias <- ifelse(df$Parameter=="cov_int_psi", ((df$X50.-Cov_alpha.psi)/Cov_alpha.psi),df$RBias)
df$RBias <- ifelse(df$Parameter=="cov_int_phi", ((df$X50.-Cov_alpha.phi)/Cov_alpha.phi),df$RBias)
df$RBias <- ifelse(df$Parameter=="cov_psi_phi", (df$X50.-Cov_psi.phi)/Cov_psi.phi, df$RBias)

#df <- filter(df, Parameter=="psi"|Parameter=="Sigma2_x"|Parameter=="Sigma2_epsilon")

# # Compute summary statistics
# summary <- df %>%
#   group_by(Q,N,ind,n_per_ind,repeats, Parameter) %>%
#   summarise(
#     mean_bias = mean(RBias),
#     .groups = 'drop'
#   )



df <- filter(df, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi"|
               Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")

df$Sim_estimate <- NA
df$Sim_estimate <- ifelse(df$Parameter=="Sigma2_intercept", Valpha, df$Sim_estimate)
df$Sim_estimate <- ifelse(df$Parameter=="Sigma2_psi",Vpsi, df$Sim_estimate)
df$Sim_estimate <- ifelse(df$Parameter=="Sigma2_phi",  Vphi, df$Sim_estimate)
df$Sim_estimate <- ifelse(df$Parameter=="cov_int_psi", Cov_alpha.psi, df$Sim_estimate)
df$Sim_estimate <- ifelse(df$Parameter=="cov_int_psi", Cov_alpha.phi, df$Sim_estimate)
df$Sim_estimate <- ifelse(df$Parameter=="cov_int_psi", Cov_psi.phi, df$Sim_estimate)

RMAD <- function(x) {
  mean_x <- mean(x)
  mae <- mean(abs(x - mean_x))
  mean_abs <- mean(abs(x))
  
  if (mean_abs == 0) return(NA_real_)
  return(mae / mean_abs)
}

# Compute summary statistics
summary <- df %>%
  group_by(Q,N,ind,n_per_ind,repeats, Parameter, Model) %>%
  summarise(
    mean_bias = mean(RBias)*100, # Single point
    sd_bias = sd(RBias),# Dispersion measure
    RMAD = RMAD(X50.)*100,
    MAE = mean(abs(RBias))*100,
    .groups = 'drop'
  )

summary_table <-  filter(summary, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi"|
                               Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")
# tables
library(xtable)
library(reshape2)


library(reshape2)

# Define parameter labels (with desired row order)
latex_labels <- c(
  "Sigma2_intercept" = "Mean $V_\\alpha$",
  "Sigma2_phi"       = "Impact $V_\\phi$",
  "Sigma2_psi"       = "Response $V_\\psi$",
  "cov_int_phi"      = "$\\text{Cov}(\\alpha, \\phi)$",
  "cov_int_psi"      = "$\\text{Cov}(\\alpha, \\psi)$",
  "cov_psi_phi"      = "$\\text{Cov}(\\psi, \\phi)$"
)
param_order <- names(latex_labels)

# Define column design order
summary_table$Design <- paste(summary_table$ind, summary_table$n_per_ind, summary_table$repeats, sep = "_")
col_order <- c("1600_2_1x","800_4_1x", "400_8_1x", "200_16_1x", "100_32_1x", "800_4_1x", "400_4_2x", "200_4_4x", "100_4_8x", "200_16_1x", "200_8_2x", "200_4_4x", "200_2_8x")

# Helper: construct LaTeX block with \multirow for metric labels
metric_block <- function(df, metric_col, label) {
  df <- df[df$Parameter %in% param_order, ]
  df$Design <- paste(df$ind, df$n_per_ind, df$repeats, sep = "_")
  
  # Aggregate (mean over replicates if needed)
  agg <- aggregate(df[[metric_col]], by = list(Parameter = df$Parameter, Design = df$Design), FUN = mean, na.rm = TRUE)
  colnames(agg)[3] <- "value"
  
  # Reshape to wide
  wide <- dcast(agg, Parameter ~ Design, value.var = "value")
  wide <- wide[match(param_order, wide$Parameter), ]
  wide$Parameter <- latex_labels[wide$Parameter]
  wide <- wide[, c("Parameter", col_order)]
  
  # Build rows with \multirow only on the first
  rows <- apply(wide, 1, function(row) {
    paste(" &", row[1], "&", paste(sprintf("%.2f", as.numeric(row[-1])), collapse = " & "), "\\\\")
  })
  rows[1] <- paste0("\\multirow{6}{*}{\\textbf{", label, "}}", rows[1])
  paste(rows, collapse = "\n")
}

# Generate each block
bias_block <- metric_block(summary_table, "mean_bias", "Bias (\\%)")
disp_block <- metric_block(summary_table, "RMAD", "Dispersion (\\%)")
mae_block  <- metric_block(summary_table, "MAE", "MAE (\\%)")

# Header rows
header1 <- " & & 1600 & 800 & 400 & 200 & 100 & 800 & 400 & 200 & 100 & 200 & 200 & 200 & 200\\\\"
header2 <- " & & 2 & 4 & 8 & 16 & 32 & 4 & 4 & 4 & 4 & 16 & 8 & 4 & 2\\\\"
header3 <- " & & 1x & 1x & 1x & 1x & 1x & 1x & 2x & 4x & 8x & 1x & 2x & 4x & 8x \\\\"

# Combine into full LaTeX table
latex_output <- paste0(
  "\\begin{table}[ht]\n",
  "\\caption{Model performance metrics (bias, dispersion, and MAE) for key variance and covariance parameters.}\n",
  "\\resizebox{\\textwidth}{!}{%\n",
  "\\begin{tabular}{llccccc|cccc|cccc}\n",
  "\\toprule\n",
  header1, "\n",
  header2, "\n",
  header3, "\n",
  "\\midrule\n",
  bias_block, "\n",
  "\\midrule\n",
  disp_block, "\n",
  "\\midrule\n",
  mae_block, "\n",
  "\\bottomrule\n",
  "\\end{tabular}%\n",
  "}\n",
  "\\end{table}\n"
)

# Print to console
cat(latex_output)


