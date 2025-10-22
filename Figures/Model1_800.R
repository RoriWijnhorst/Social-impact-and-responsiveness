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

df2_800_400_2_1x  <-read.csv("Output/res_800_400_2_1x_M2.csv")
df2_800_200_4_1x  <-read.csv("Output/res_800_200_4_1x_M2.csv")
df2_800_100_8_1x  <-read.csv("Output/res_800_100_8_1x_M2.csv")
df2_800_50_16_1x  <-read.csv("Output/res_800_50_16_1x_M2.csv")

df3_800_400_2_1x  <-read.csv("Output/res_800_400_2_1x_M2.csv")
df3_800_200_4_1x  <-read.csv("Output/res_800_200_4_1x_M2.csv")
df3_800_100_8_1x  <-read.csv("Output/res_800_100_8_1x_M2.csv")
df3_800_50_16_1x  <-read.csv("Output/res_800_50_16_1x_M2.csv")


df <- gdata::combine(
  df1_800_400_2_1x, 
  df1_800_200_4_1x, 
  df1_800_100_8_1x,
  df1_800_50_16_1x,
  
  df2_800_400_2_1x, 
  df2_800_200_4_1x, 
  df2_800_100_8_1x,
  df2_800_50_16_1x,
  
  df3_800_400_2_1x, 
  df3_800_200_4_1x, 
  df3_800_100_8_1x,
  df3_800_50_16_1x
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
df <- df %>%  mutate(ind = fct_relevel(ind, "400","200","100","50"))

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
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Bias (%)", limits=c(-15,15), breaks = seq(-15, 15, 5)) +
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


p2<-ggplot(summary1, aes(x = as.numeric(ind), y = RMAD, color = Parameter,  group = Model)) + 
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=RMAD,linetype=Model), size=0.7) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Dispersion (%)", limits=c(10,40), breaks = seq(10, 40, 5)) +
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
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Bias (%)", limits=c(-15,15), breaks = seq(-15, 15, 5)) +
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
        legend.position = "none",
        legend.text=element_text(size=10),
        legend.title = element_text())
  

p4<-ggplot(summary2, aes(x = as.numeric(ind), y = RMAD, color = Parameter,  group = Model)) + 
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=RMAD,linetype=Model), size=0.7) +
  scale_x_continuous(name="number of individuals\nsocial partners", 
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Dispersion (%)", limits=c(10,40), breaks = seq(10, 40, 5)) +
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
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
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
                     breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16")) +
  scale_y_continuous(name="Bias (%)", limits=c(10,45), breaks = seq(10, 45, 5)) +
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


## Create plot ----
patch1 <- ((p1 + p3 + p5) / (p2 + p4 + p6)) + plot_layout(guides = 'collect')
wrap_elements(panel = patch1) + 
  labs(tag = "Number of individuals/\nsocial partners") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  ) 

