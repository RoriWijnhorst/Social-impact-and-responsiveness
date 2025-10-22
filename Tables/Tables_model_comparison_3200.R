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

df2_3200_1600_2_1x  <-read.csv("Output/res_3200_1600_2_1x_M2.csv") 
df2_3200_800_4_1x   <-read.csv("Output/res_3200_800_4_1x_M2.csv") 
df2_3200_400_8_1x   <-read.csv("Output/res_3200_400_8_1x_M2.csv") 
df2_3200_200_16_1x  <-read.csv("Output/res_3200_200_16_1x_M2.csv")
df2_3200_100_32_1x  <-read.csv("Output/res_3200_100_32_1x_M2.csv")

df3_3200_1600_2_1x  <-read.csv("Output/res_3200_1600_2_1x_M3.csv") 
df3_3200_800_4_1x   <-read.csv("Output/res_3200_800_4_1x_M3.csv") 
df3_3200_400_8_1x   <-read.csv("Output/res_3200_400_8_1x_M3.csv") 
df3_3200_200_16_1x  <-read.csv("Output/res_3200_200_16_1x_M3.csv")
df3_3200_100_32_1x  <-read.csv("Output/res_3200_100_32_1x_M3.csv")


df4_3200_1600_2_1x  <-read.csv("Output/res_3200_1600_2_1x_M4.csv") 
df4_3200_800_4_1x   <-read.csv("Output/res_3200_800_4_1x_M4.csv") 
df4_3200_400_8_1x   <-read.csv("Output/res_3200_400_8_1x_M4.csv") 
df4_3200_200_16_1x  <-read.csv("Output/res_3200_200_16_1x_M4.csv")
df4_3200_100_32_1x  <-read.csv("Output/res_3200_100_32_1x_M4.csv")

df5_3200_1600_2_1x  <-read.csv("Output/res_3200_1600_2_1x_M5.csv") 
df5_3200_800_4_1x   <-read.csv("Output/res_3200_800_4_1x_M5.csv") 
df5_3200_400_8_1x   <-read.csv("Output/res_3200_400_8_1x_M5.csv") 
df5_3200_200_16_1x  <-read.csv("Output/res_3200_200_16_1x_M5.csv")
df5_3200_100_32_1x  <-read.csv("Output/res_3200_100_32_1x_M5.csv")

df <- gdata::combine(
  df1_3200_1600_2_1x, 
  df1_3200_800_4_1x, 
  df1_3200_400_8_1x, 
  df1_3200_200_16_1x,
  df1_3200_100_32_1x,
  
  df2_3200_1600_2_1x, 
  df2_3200_800_4_1x, 
  df2_3200_400_8_1x, 
  df2_3200_200_16_1x,
  df2_3200_100_32_1x,
  
  df3_3200_1600_2_1x, 
  df3_3200_800_4_1x, 
  df3_3200_400_8_1x, 
  df3_3200_200_16_1x,
  df3_3200_100_32_1x,
  
  df4_3200_1600_2_1x, 
  df4_3200_800_4_1x, 
  df4_3200_400_8_1x, 
  df4_3200_200_16_1x,
  df4_3200_100_32_1x,
  
  df5_3200_1600_2_1x, 
  df5_3200_800_4_1x, 
  df5_3200_400_8_1x, 
  df5_3200_200_16_1x,
  df5_3200_100_32_1x
  )

# Means
Mpsi     = 0.3    # interaction coefficient psi

# Variances
Valpha   = 0.2    # variance direct effect (mean behaviour)
Vepsilon = 0.01   # variance social partner effect (residual impact)
Vx       = 1      # variance social partner phenotype (impact covariate)

# Correlations between variance components
r_alpha_epsilon  =  0       # part of Cov int-phi
r_alpha_x        =  0.6     # part of Cov int-phi 

Vphi <- Mpsi^2 * Vx + Vepsilon
Cov_alpha.phi <- r_alpha_epsilon*(sqrt(Valpha*Vepsilon))+ Mpsi*r_alpha_x*(sqrt(Valpha*Vx))


df$RBias <- NA
df$RBias <- ifelse(df$Parameter=="psi", (df$X50.-Mpsi)/Mpsi, df$RBias)
df$RBias <- ifelse(df$Parameter=="Sigma2_phi", ((df$X50.-Vphi)/Vphi),df$RBias)
df$RBias <- ifelse(df$Parameter=="cov_int_phi", ((df$X50.-Cov_alpha.phi)/Cov_alpha.phi),df$RBias)

df <- rename(df, Q=source)
df$Q <- gsub('df','',df$Q)
df$Q <- as.factor(df$Q)
df[c('Model','N','ind','n_per_ind','repeats')] <- str_split_fixed(df$Q, '_', 5)
df <- df %>%  mutate(ind = fct_relevel(ind, "1600", "800", "400","200","100"))

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

summary_table_psi <-  filter(summary, Parameter ==  "psi")
summary_table_var <-  filter(summary, Parameter ==  "Sigma2_phi")
summary_table_cov <-  filter(summary, Parameter ==  "cov_int_phi")


# tables
library(xtable)
library(reshape2)

df<-summary_table_var

# Reformat factor levels
df$Model <- factor(df$Model, labels = c("I\\&R", "V--P", "Trait", "Alt1", "Alt2"))
df$ind <- factor(df$ind, levels = c(1600, 800, 400, 200, 100), ordered = TRUE)

# Function to reshape and round
reshape_metric <- function(df, metric_name) {
  dcast(df, ind ~ Model, value.var = metric_name)
}

# Generate tables
bias <- reshape_metric(df, "mean_bias")
disp <- reshape_metric(df, "RMAD")

# Combine all with labels
format_block <- function(metric_df, label) {
  apply(metric_df, 1, function(row) {
    sprintf("   & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\",
            row[1], as.numeric(row[2]), as.numeric(row[3]), as.numeric(row[4]),
            as.numeric(row[5]), as.numeric(row[6]))
  }) |> paste(collapse = "\n")
}


# Compose final LaTeX
latex_output <- paste0(
  "\\begin{table}[ht]\n",
  "\\begin{tabular}{lcccccc}\n",
  "  \\toprule \n",
  " & Individuals & I\\&R & V--P & Trait & Trait+EIV &Trait+RS\\\\ \n",
  "  \\midrule\n",
  "\\multirow{5}{*}{\\textbf{Bias (\\%)}} \n", format_block(bias, "Bias (%)"), "\n",
  "  \\midrule\n",
  "\\multirow{5}{*}{\\textbf{Dispersion (\\%)}} \n", format_block(disp, "Dispersion (%)"), "\n",
  "  \\midrule\n",
  "  \\end{tabular}\n",
  "\\end{table}"
)
# Print it
cat(latex_output)


df<-summary_table_cov

# Reformat factor levels
df$Model <- factor(df$Model, labels = c("I\\&R", "V--P", "Trait", "Alt1", "Alt2"))
df$ind <- factor(df$ind, levels = c(1600, 800, 400, 200, 100), ordered = TRUE)

# Function to reshape and round
reshape_metric <- function(df, metric_name) {
  dcast(df, ind ~ Model, value.var = metric_name)
}

# Generate tables
bias <- reshape_metric(df, "mean_bias")
disp <- reshape_metric(df, "RMAD")

# Combine all with labels
format_block <- function(metric_df, label) {
  apply(metric_df, 1, function(row) {
    sprintf("   & %s & %.2f & %.2f & %.2f & %.2f & %.2f \\\\",
            row[1], as.numeric(row[2]), as.numeric(row[3]), as.numeric(row[4]),
            as.numeric(row[5]), as.numeric(row[6]))
  }) |> paste(collapse = "\n")
}


# Compose final LaTeX
latex_output <- paste0(
  "\\begin{table}[ht]\n",
  "\\begin{tabular}{lcccccc}\n",
  "  \\toprule \n",
  " & Individuals & I\\&R & V--P & Trait & Trait+EIV &Trait+RS\\\\ \n",
  "  \\midrule\n",
  "\\multirow{5}{*}{\\textbf{Bias (\\%)}} \n", format_block(bias, "Bias (%)"), "\n",
  "  \\midrule\n",
  "\\multirow{5}{*}{\\textbf{Dispersion (\\%)}} \n", format_block(disp, "Dispersion (%)"), "\n",
  "  \\midrule\n",
  "  \\end{tabular}\n",
  "\\end{table}"
)
# Print it
cat(latex_output)



## table format
df<-summary_table_psi

# Reformat factor levels
df$Model <- factor(df$Model)
df$ind <- factor(df$ind, levels = c(1600, 800, 400, 200, 100), ordered = TRUE)

# Function to reshape and round
reshape_metric <- function(df, metric_name) {
  dcast(df, ind ~ Model, value.var = metric_name)
}

# Generate tables
bias <- reshape_metric(df, "mean_bias")
disp <- reshape_metric(df, "RMAD")

# Combine all with labels
# Combine all with labels
format_block <- function(metric_df, label) {
  apply(metric_df, 1, function(row) {
    sprintf("   & %s & %.2f & %.2f & %.2f & %.2f \\\\",
            row[1], as.numeric(row[2]), as.numeric(row[3]), as.numeric(row[4]),
            as.numeric(row[5]))
  }) |> paste(collapse = "\n")
}

# Compose final LaTeX
latex_output <- paste0(
  "\\begin{table}[ht]\n",
  "\\begin{tabular}{lccccc}\n",
  "  \\toprule \n",
  " & Individuals & I\\&R & Trait & Trait+EIV & Trait+RS\\\\ \n",
  "  \\midrule\n",
  "\\multirow{5}{*}{\\textbf{Bias (\\%)}} \n", format_block(bias, "Bias (%)"), "\n",
  "  \\midrule\n",
  "\\multirow{5}{*}{\\textbf{Dispersion (\\%)}} \n", format_block(disp, "Dispersion (%)"), "\n",
  "  \\midrule\n",
  "  \\end{tabular}\n",
  "\\end{table}"
)

# Print it
cat(latex_output)

