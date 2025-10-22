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

df <- gdata::combine(
  df1_3200_1600_2_1x, 
  df1_3200_800_4_1x ,
  df1_3200_400_8_1x ,
  df1_3200_200_16_1x,
  df1_3200_100_32_1x,
  
  df2_3200_1600_2_1x, 
  df2_3200_800_4_1x ,
  df2_3200_400_8_1x ,
  df2_3200_200_16_1x,
  df2_3200_100_32_1x,
  
  df3_3200_1600_2_1x, 
  df3_3200_800_4_1x ,
  df3_3200_400_8_1x ,
  df3_3200_200_16_1x,
  df3_3200_100_32_1x
  )

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
df[c('Model','N','ind','n_per_ind','repeats')] <- str_split_fixed(df$Q, '_', 6)
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



RMAD <- function(x) {
  mean_x <- mean(x)
  mae <- mean(abs(x - mean_x))
  mean_abs <- mean(abs(x))
  
  if (mean_abs == 0) return(NA_real_)
  return(mae / mean_abs)
}

# Compute summary statistics
summary <- df %>%
  group_by(Q,N,ind,n_per_ind,repeats, Parameter,Model) %>%
  summarise(
    mean_bias = mean(RBias)*100, # Single point
    sd_bias = sd(RBias)*100,# Dispersion measure
    RMAD = RMAD(X50.)*100,
    .groups = 'drop'
  )


summary_table_psi <-  filter(summary, Parameter ==  "psi")
summary_table_var <-  filter(summary, Parameter ==  "Sigma2_phi")
summary_table_cov <-  filter(summary, Parameter ==  "cov_int_phi")

summary1 <- filter(summary, Parameter == "psi")
summary2 <- filter(summary, Parameter == "Sigma2_phi" )
summary3 <- filter(summary, Parameter == "cov_int_phi")

summary1$ind <- as.numeric(as.factor(summary1$ind))

# === Bias Plot ===
p1 <- ggplot(summary1, aes(x = as.numeric(ind), y = mean_bias, color = Parameter,  group = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size=0.4) +  # Reference line at 0
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=mean_bias,linetype=Model), size=0.7) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:5), labels=c("1600\n2","800\n4","400\n8","200\n16", "100\n32")) +
  scale_y_continuous(name="Bias (%)", limits=c(-20,10), breaks = seq(-20, 10, 5)) +
  scale_colour_manual(values=c("#E5C494"), labels=c("\U1D449\U1D719"),
                      guide="none") +
  scale_shape_discrete(labels=c("I&R", "V-P", "Trait"))+
  scale_linetype_discrete(labels=c("I&R", "V-P", "Trait"))+
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(),         
        axis.text = element_text(color="black"),
        axis.title.x=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=10),
        legend.title = element_text())


p2<-ggplot(summary1, aes(x = as.numeric(ind), y = RMAD, color = Parameter,  group = Model)) + 
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=RMAD,linetype=Model), size=0.7) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:5), labels=c("1600\n2","800\n4","400\n8","200\n16", "100\n32")) +
  scale_y_continuous(name="Dispersion (%)", limits=c(0,26), breaks = seq(0, 25, 5)) +
  scale_colour_manual(values=c("#E5C494"), labels=c("\U1D449\U1D719"),
                      guide="none") +
  scale_shape_discrete(labels=c("I&R", "V-P", "Trait"))+
  scale_linetype_discrete(labels=c("I&R", "V-P", "Trait"))+
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),         
        axis.text = element_text(color="black"),
        axis.title.x=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=10),
        legend.title = element_text())


p3 <- ggplot(summary2, aes(x = as.numeric(ind), y = mean_bias, color = Parameter,  group = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size=0.4) +  # Reference line at 0
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=mean_bias,linetype=Model), size=0.7) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:5), labels=c("1600\n2","800\n4","400\n8","200\n16", "100\n32")) +
  scale_y_continuous(name="Bias (%)", limits=c(-20,10), breaks = seq(-20, 10, 5)) +
  scale_colour_manual(values=c("#66C2A5"), labels=c("\U1D449\U1D719"),
                      guide="none") +
  scale_shape_discrete(labels=c("I&R", "V-P", "Trait"))+
  scale_linetype_discrete(labels=c("I&R", "V-P", "Trait"))+
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),         
        axis.text.x = element_blank(),
        axis.text = element_text(color="black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=10),
        legend.title = element_text())
  

p4<-ggplot(summary2, aes(x = as.numeric(ind), y = RMAD, color = Parameter,  group = Model)) + 
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=RMAD,linetype=Model), size=0.7) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:5), labels=c("1600\n2","800\n4","400\n8","200\n16", "100\n32")) +
  scale_y_continuous(name="Dispersion (%)", limits=c(0,26), breaks = seq(0, 25, 5)) +
  scale_colour_manual(values=c("#66C2A5"), labels=c("\U1D449\U1D719"),
                      guide="none") +
  scale_shape_discrete(labels=c("I&R", "V-P", "Trait"))+
  scale_linetype_discrete(labels=c("I&R", "V-P", "Trait"))+
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),         
        axis.text = element_text(color="black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "right",
        legend.text=element_text(size=10),
        legend.title = element_text())


p5 <- ggplot(summary3, aes(x = as.numeric(ind), y = mean_bias, color = Parameter,  group = Model)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size=0.4) +  # Reference line at 0
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=mean_bias,linetype=Model), size=0.7) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:5), labels=c("1600\n2","800\n4","400\n8","200\n16", "100\n32")) +
  scale_y_continuous(name="Bias (%)", limits=c(-65,2), breaks = seq(-65, 2, 5)) +
  scale_colour_manual(values=c("#3F4176"), labels=c("\U1D449\U1D719"),
                      guide="none") +
  scale_shape_discrete(labels=c("I&R", "V-P", "Trait"))+
  scale_linetype_discrete(labels=c("I&R", "V-P", "Trait"))+
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),         
        axis.text.x = element_blank(),
        axis.text = element_text(color="black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "right",
        legend.text=element_text(size=10),
        legend.title = element_text())


p6<-ggplot(summary3, aes(x = as.numeric(ind), y = RMAD, color = Parameter,  group = Model)) + 
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=RMAD,linetype=Model), size=0.7) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:5), labels=c("1600\n2","800\n4","400\n8","200\n16", "100\n32")) +
  scale_y_continuous(name="Dispersion (%)", limits=c(0,26), breaks = seq(0, 25, 5)) +
  scale_colour_manual(values=c("#3F4176"), labels=c("\U1D449\U1D719"),
                      guide="none") +
  scale_shape_discrete(labels=c("I&R", "V-P", "Trait"))+
  scale_linetype_discrete(labels=c("I&R", "V-P", "Trait"))+
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),         
        axis.text = element_text(color="black"),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "right",
        legend.text=element_text(size=10),
        legend.title = element_text())


# Three title "plots"
t1 <- ggplot() + theme_void() +
  labs(title = ("Population response \U1D713")) +
  theme(plot.title = element_text(hjust = 0.4, size = 12, margin = margin(b = 2)))

t2 <- ggplot() + theme_void() +
  labs(title = ("Social impact \U1D449\U1D719")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12, margin = margin(b = 2)))

t3 <- ggplot() + theme_void() +
  labs(title = ("Mean - impact \U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D719\U2769")) +
  theme(plot.title = element_text(hjust = 0.7, size = 12, margin = margin(b = 2)))

top_titles <- t1 + t2 + t3  # a 1×3 row

# Your existing panels
main <- ((p1 + p3 + p5) / (p2 + p4 + p6)) + plot_layout(guides = "collect")

xlab_row <- ggplot() + theme_void() +
  labs(title = "Number of individuals\nsocial partners") +
  theme(plot.title = element_text(hjust = 0.5, size = 11,
                                  margin = margin(t = 2, b = 4)))

final_plot <- (top_titles / main / xlab_row) +
  plot_layout(heights = c(0, 1, 0.01))
# export pdf 9 x5.5 inch
final_plot

diff_to_zero <- function(x) {
  ref <- which.min(abs(x))
  round(x - x[ref], 2)
}

df<-summary_table_var

# Setup
ind_vals <- sort(unique(df$ind), decreasing = TRUE)  # 400 to 50
models <- paste0("Model ", 1:3)

# Initialize
bias_mat <- matrix(NA, nrow = 3, ncol = length(ind_vals))
rmad_mat <- matrix(NA, nrow = 3, ncol = length(ind_vals))
mae_mat <- matrix(NA, nrow = 3, ncol = length(ind_vals))

colnames(bias_mat) <- colnames(rmad_mat) <- colnames(mae_mat) <- as.character(ind_vals)
rownames(bias_mat) <- rownames(rmad_mat) <- rownames(mae_mat) <- models

# Fill in matrices
for (i in seq_along(ind_vals)) {
  ind <- ind_vals[i]
  sub <- df[df$ind == ind, ]
  bias_mat[, i] <- diff_to_zero(sub$mean_bias)
  rmad_mat[, i] <- diff_to_zero(sub$RMAD)
  mae_mat[, i] <- diff_to_zero(sub$MAE)
}

# Convert to data frames and label
bias_df <- data.frame(Model = models, Metric = "Bias (%)", bias_mat)
rmad_df <- data.frame(Model = models, Metric = "Dispersion (%)", rmad_mat)
mae_df  <- data.frame(Model = models, Metric = "MAE (%)", mae_mat)

# Combine all
summary_df <- rbind(bias_df, rmad_df, mae_df)

df<-summary_table_var

# Setup
ind_vals <- sort(unique(df$ind), decreasing = TRUE)  # 400 to 50
models <- paste0("Model ", 1:3)

# Initialize
bias_mat <- matrix(NA, nrow = 3, ncol = length(ind_vals))
rmad_mat <- matrix(NA, nrow = 3, ncol = length(ind_vals))
mae_mat <- matrix(NA, nrow = 3, ncol = length(ind_vals))

colnames(bias_mat) <- colnames(rmad_mat) <- colnames(mae_mat) <- as.character(ind_vals)
rownames(bias_mat) <- rownames(rmad_mat) <- rownames(mae_mat) <- models

# Fill in matrices
for (i in seq_along(ind_vals)) {
  ind <- ind_vals[i]
  sub <- df[df$ind == ind, ]
  bias_mat[, i] <- diff_to_zero(sub$mean_bias)
  rmad_mat[, i] <- diff_to_zero(sub$RMAD)
  mae_mat[, i] <- diff_to_zero(sub$MAE)
}

# Convert to data frames and label
bias_df <- data.frame(Model = models, Metric = "Bias (%)", bias_mat)
rmad_df <- data.frame(Model = models, Metric = "Dispersion (%)", rmad_mat)
mae_df  <- data.frame(Model = models, Metric = "MAE (%)", mae_mat)

# Combine all
summary_df <- rbind(bias_df, rmad_df, mae_df)


# Save to CSV
write.csv(summary_df, "model_comparison_summary.csv", row.names = FALSE)


install.packages("xtable")  # if not already installed
library(xtable)
df<-summary_table_cov

# Reformat factor levels
df$Model <- factor(df$Model, labels = c("I\\&R", "V--P", "Trait"))
df$ind <- factor(df$ind, levels = c(400, 200, 100, 50), ordered = TRUE)

# Function to reshape and round
reshape_metric <- function(df, metric_name) {
  dcast(df, ind ~ Model, value.var = metric_name)
}

# Generate tables
bias <- reshape_metric(df, "mean_bias")
disp <- reshape_metric(df, "RMAD")
mae  <- reshape_metric(df, "MAE")

# Combine all with labels
format_block <- function(metric_df, label) {
  apply(metric_df, 1, function(row) {
    sprintf("   & %s & %.2f & %.2f & %.2f \\\\", row[1], as.numeric(row[2]), as.numeric(row[3]), as.numeric(row[4]))
  }) |> paste(collapse = "\n")
}

# Compose final LaTeX
latex_output <- paste0(
  "\\begin{table}[ht]\n",
  "\\begin{tabular}{lcccc}\n",
  "\\multicolumn{5}{c}{\\textbf{Covariance: mean behaviour-social impact}}\\\\\n",
  "  \\toprule \n",
  " & Individuals & I\\&R & V--P & Trait \\\\ \n",
  "  \\midrule\n",
  "\\multirow{4}{*}{\\textbf{Bias (\\%)}} \n", format_block(bias, "Bias (%)"), "\n",
  "  \\midrule\n",
  "\\multirow{4}{*}{\\textbf{Dispersion (\\%)}} \n", format_block(disp, "Dispersion (%)"), "\n",
  "  \\midrule\n",
  "\\multirow{4}{*}{\\textbf{MAE (\\%)}} \n", format_block(mae, "MAE (%)"), "\n",
  "  \\bottomrule\n",
  "  \\end{tabular}\n",
  "\\end{table}"
)

# Print it
cat(latex_output)



## table format
df<-summary_table_psi

# Reformat factor levels
df$Model <- factor(df$Model, labels = c("I\\&R", "Trait"))
df$ind <- factor(df$ind, levels = c(400, 200, 100, 50), ordered = TRUE)

# Function to reshape and round
reshape_metric <- function(df, metric_name) {
  dcast(df, ind ~ Model, value.var = metric_name)
}

# Generate tables
bias <- reshape_metric(df, "mean_bias")
disp <- reshape_metric(df, "RMAD")
mae  <- reshape_metric(df, "MAE")

# Combine all with labels
format_block <- function(metric_df, label) {
  apply(metric_df, 1, function(row) {
    sprintf("   & %s & %.2f & %.2f \\\\", row[1], as.numeric(row[2]), as.numeric(row[3]))
  }) |> paste(collapse = "\n")
}

# Compose final LaTeX
latex_output <- paste0(
  "\\begin{table}[ht]\n",
  "\\begin{tabular}{lccc}\n",
  "\\multicolumn{4}{c}{\\textbf{Covariance: mean behaviour-social impact}}\\\\\n",
  "  \\toprule \n",
  " & Individuals & I\\&R & Trait\\\\ \n",
  "  \\midrule\n",
  "\\multirow{4}{*}{\\textbf{Bias (\\%)}} \n", format_block(bias, "Bias (%)"), "\n",
  "  \\midrule\n",
  "\\multirow{4}{*}{\\textbf{Dispersion (\\%)}} \n", format_block(disp, "Dispersion (%)"), "\n",
  "  \\midrule\n",
  "\\multirow{4}{*}{\\textbf{MAE (\\%)}} \n", format_block(mae, "MAE (%)"), "\n",
  "  \\bottomrule\n",
  "  \\end{tabular}\n",
  "\\end{table}"
)

# Print it
cat(latex_output)