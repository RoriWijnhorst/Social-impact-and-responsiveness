library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(forcats)


df1_800_400_2_1x  <-read.csv("Output/res_800_400_2_1x_M1.csv") 
df1_800_200_4_1x  <-read.csv("Output/res_800_200_4_1x_M1.csv")
df1_800_100_8_1x  <-read.csv("Output/res_800_100_8_1x_M1.csv")
df1_800_50_16_1x  <-read.csv("Output/res_800_50_16_1x_M1.csv")

df2_800_200_4_1x <-read.csv("Output/res_800_200_4_1x_M1.csv")
df2_800_100_4_2x <-read.csv("Output/res_800_100_4_2x_M1.csv")
df2_800_50_4_4x  <-read.csv("Output/res_800_50_4_4x_M1.csv")

df3_800_100_8_1x <-read.csv("Output/res_800_100_8_1x_M1.csv")
df3_800_100_4_2x <-read.csv("Output/res_800_100_4_2x_M1.csv")
df3_800_100_2_4x <-read.csv("Output/res_800_100_2_4x_M1.csv")


df <- gdata::combine(
  df1_800_400_2_1x,
  df1_800_200_4_1x,
  df1_800_100_8_1x,
  df1_800_50_16_1x,
  
  df2_800_200_4_1x,
  df2_800_100_4_2x,
  df2_800_50_4_4x,
  
  df3_800_100_8_1x,
  df3_800_100_4_2x,
  df3_800_100_2_4x)

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
df <- df %>%  mutate(ind = fct_relevel(ind, "400","200","100","50"))

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
  group_by(Q,N,ind,n_per_ind,repeats, Parameter, Sim_estimate) %>%
  summarise(
    mean_bias = mean(RBias), # Single point
    sd_bias = sd(RBias),# Dispersion measure
    RMAD = RMAD(X50.),
    .groups = 'drop'
  )

# Plot bias and dispersion
dfv_Q1 <- filter(summary, Q==1, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi")
dfv_Q1$ind <- as.numeric(as.factor(dfv_Q1$ind))

dfv_Q2 <- filter(summary, Q==2, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi")
dfv_Q2$ind <- as.numeric(as.factor(dfv_Q2$ind))


dfv_Q3 <- filter(summary, Q==3, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi")
dfv_Q3$n_per_ind <- fct_relevel(dfv_Q3$n_per_ind, "8","4","2")
dfv_Q3$n_per_ind <- as.numeric(as.factor(dfv_Q3$n_per_ind))

# Plot bias and dispersion
dfc_Q1 <- filter(summary, Q==1, Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")
dfc_Q1$ind <- as.numeric(as.factor(dfc_Q1$ind))

dfc_Q2 <- filter(summary, Q==2, Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")
dfc_Q2$ind <- as.numeric(as.factor(dfc_Q2$ind))

dfc_Q3 <- filter(summary, Q==3, Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")
dfc_Q3$n_per_ind <- fct_relevel(dfc_Q3$n_per_ind, "8","4","2")
dfc_Q3$n_per_ind <- as.numeric(as.factor(dfc_Q3$n_per_ind))



# === Bias Plot ===
B1 <- ggplot(dfv_Q1, aes(x = as.numeric(ind), y = mean_bias*100, color = Parameter, , group = Parameter)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size=0.4) +  # Reference line at 0
  geom_line(size=0.4) + 
  geom_point(size=1) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Bias (%)",limits=c(-5,15), breaks = seq(-5, 15, 5)) +
  scale_color_manual(values=c("#2E6F5D", "#66C2A5", "#B2E1D3"), 
                     labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
  annotate("text", x=2.5, y=15, label= "repeats = 1",size=3) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_line(),
    axis.title.x = element_blank(),
    axis.text = element_text(color="black"),
    legend.position = "none"
  ) 

D1 <- ggplot(dfv_Q1, aes(x = ind, y = RMAD*100, color = Parameter, group = Parameter)) + 
  geom_line(size=0.4) +
  geom_point(size=1) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Dispersion (%)", limits=c(15, 37)) +
  scale_color_manual(values=c("#2E6F5D", "#66C2A5", "#B2E1D3")) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(color="black"),
    legend.position = "none"
  ) 

# === Bias Plot ===
B2 <- ggplot(dfv_Q2, aes(x = as.numeric(ind), y = mean_bias*100, color = Parameter, group = Parameter)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size=0.4) + # Reference line at 0
  geom_line(size=0.4) + 
  geom_point(size=1) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Bias (%)", limits=c(-5,15), breaks = seq(-5, 15, 5)) +
  scale_color_manual(values=c("#2E6F5D", "#66C2A5", "#B2E1D3"), 
                     labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
  theme_bw(base_size = 8) +
  annotate("text", x=3, y=15, label= "partners = 4",size=3) +
  labs(title="             Variances")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(face = "bold"),
    axis.text = element_text(color="black"),
    legend.position = "none"
  ) 

D2 <- ggplot(dfv_Q2, aes(x = ind, y = RMAD*100, color = Parameter, group = Parameter)) + 
  geom_line(size=0.4) +
  geom_point(size=1) +
  scale_x_continuous(name="number of individuals\nrepeats", breaks=c(2,3,4), labels=c("200\n1x","100\n2x","50\n4x"))+
  scale_y_continuous(name="Dispersion (%)", limits=c(15, 37)) +
  scale_color_manual(values=c("#2E6F5D", "#66C2A5", "#B2E1D3")) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.y = element_blank(),
    axis.text = element_text(color="black"),
    axis.ticks = element_line(),
    legend.position = "none"
  ) 


# === Bias Plot ===
B3 <- ggplot(dfv_Q3, aes(x = as.numeric(n_per_ind), y = mean_bias*100, color = Parameter, group = Parameter)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size=0.4) +  # Reference line at 0
  geom_line(size=0.4) + 
  geom_point(size=1) +
  scale_x_continuous(name="number of social partners\nrepeats", breaks=c(1,2,3), labels=c("8\n1x","4\n2x","2\n4x"))+
  scale_y_continuous(name="Bias (%)",limits=c(-5,15), breaks = seq(-5, 15, 5)) +
  scale_color_manual(values=c("#2E6F5D", "#66C2A5", "#B2E1D3"), 
                     labels=c("V_int", "V_phi", "V_psi")) +
  annotate("text", x=2, y=15, label= "individuals = 100",size=3) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(color="black"),
    legend.position = "right"
  ) 

D3 <- ggplot(dfv_Q3, aes(x = as.numeric(n_per_ind), y = RMAD*100, color = Parameter, group = Parameter)) + 
  geom_line(size=0.4) + 
  geom_point(size=1) +
  scale_x_continuous(name="number of social partners\nrepeats", breaks=c(1,2,3), labels=c("8\n1x","4\n2x","2\n4x"))+
  scale_y_continuous(name="Dispersion (%)", limits=c(15, 37)) +
  scale_color_manual(values=c("#2E6F5D", "#66C2A5", "#B2E1D3")) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.y = element_blank(),
    axis.text = element_text(color="black"),
    axis.ticks = element_line(),
    legend.position = "none"
  ) 



# === Bias Plot ===
C_b1 <- ggplot(dfc_Q1, aes(x = as.numeric(ind), y = mean_bias*100, color = Parameter, , group = Parameter)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size=0.4) +  # Reference line at 0
  geom_line(size=0.4) + 
  geom_point(size=1) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Bias (%)", limits=c(-35,5), breaks = seq(-30,0,10))+
  scale_color_manual(values=c( "#3F4176", "#6D6EB0", "#C0C0DC"), 
                     labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
  annotate("text", x=2.5, y=4.4, label= "repeats = 1",size=3) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_line(),
    axis.title.x = element_blank(),
    axis.text = element_text(color="black"),
    legend.position = "none"
  ) 

C_d1 <- ggplot(dfc_Q1, aes(x = ind, y = RMAD*100, color = Parameter, group = Parameter)) + 
  geom_line(size=0.4) +
  geom_point(size=1) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Dispersion (%)", limits=c(15, 37)) +
  scale_color_manual(values=c( "#3F4176", "#6D6EB0", "#C0C0DC")) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(color="black"),
    legend.position = "none"
  ) 

# Plot bias and dispersion


# === Bias Plot ===
C_b2 <- ggplot(dfc_Q2, aes(x = as.numeric(ind), y = mean_bias*100, color = Parameter, group = Parameter)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size=0.4) +  # Reference line at 0
  geom_line(size=0.4) + 
  geom_point(size=1) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Bias (%)", limits=c(-35,5), breaks = seq(-30,0,10))+
  scale_color_manual(values=c( "#3F4176", "#6D6EB0", "#C0C0DC"), 
                     labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
  theme_bw(base_size = 8) +
  annotate("text", x=3, y=4.4, label= "partners = 4",size=3) +
  labs(title="             Covariances")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(face = "bold"),
    axis.text = element_text(color="black"),
    legend.position = "none"
  ) 

C_d2 <- ggplot(dfc_Q2, aes(x = ind, y = RMAD*100, color = Parameter, group = Parameter)) + 
  geom_line(size=0.4) +
  geom_point(size=1) +
  scale_x_continuous(name="number of individuals\nrepeats", breaks=c(2,3,4), labels=c("200\n1x","100\n2x","50\n4x"))+
  scale_y_continuous(name="Dispersion (%)", limits=c(15, 37)) +
  scale_color_manual(values=c( "#3F4176", "#6D6EB0", "#C0C0DC")) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.y = element_blank(),
    axis.text = element_text(color="black"),
    axis.ticks = element_line(),
    legend.position = "none"
  ) 

# Plot bias and dispersion

# === Bias Plot ===
C_b3 <- ggplot(dfc_Q3, aes(x = as.numeric(n_per_ind), y = mean_bias*100, color = Parameter, group = Parameter)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", size=0.4) +  # Reference line at 0
  geom_line(size=0.4) + 
  geom_point(size=1) +
  scale_x_continuous(name="number of social partners\nrepeats", breaks=c(1,2,3), labels=c("8\n1x","4\n2x","2\n4x"))+
  scale_y_continuous(name="Bias (%)",limits=c(-35,5), breaks = seq(-30,0,10)) +
  scale_color_manual(values=c( "#3F4176", "#6D6EB0", "#C0C0DC"), 
                     labels=c("Cov_int_phi", "Cov_int_psi", "Cov_psi_phi")) +
  
  annotate("text", x=2, y=4.4, label= "individuals = 100",size=3) +
  
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(color="black"),
    legend.position = "right"
  ) 

C_d3 <- ggplot(dfc_Q3, aes(x = as.numeric(n_per_ind), y = RMAD*100, color = Parameter, group = Parameter)) + 
  geom_line(size=0.4) + 
  geom_point(size=1) +
  scale_x_continuous(name="number of social partners\nrepeats", breaks=c(1,2,3), labels=c("8\n1x","4\n2x","2\n4x"))+
  scale_y_continuous(name="Dispersion (%)", limits=c(15, 37)) +
  scale_color_manual(values=c( "#3F4176", "#6D6EB0", "#C0C0DC")) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.y = element_blank(),
    axis.text = element_text(color="black"),
    axis.ticks = element_line(),
    legend.position = "none"
  ) 

patch1 <- (B1+B2+B3)/(D1 + D2 + D3)/(C_b1+C_b2+C_b3)/(C_d1 + C_d2 + C_d3) + plot_layout(guides = 'collect')
wrap_elements(panel = patch1) + 
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  ) 

