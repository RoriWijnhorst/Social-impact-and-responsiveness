library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(forcats)

setwd("~/Social-impact-and-responsiveness/Model_comparison")

df1_800_400_2_1x  <-read.csv("res_800_400_2_1x_M1.csv")
df1_800_200_4_1x  <-read.csv("res_800_200_4_1x_M1.csv")
df1_800_100_8_1x  <-read.csv("res_800_100_8_1x_M1.csv")
df1_800_50_16_1x  <-read.csv("res_800_50_16_1x_M1.csv")

df4_800_400_2_1x  <-read.csv("res_800_400_2_1x_M4.csv")
df4_800_200_4_1x  <-read.csv("res_800_200_4_1x_M4.csv")
df4_800_100_8_1x  <-read.csv("res_800_100_8_1x_M4.csv")
df4_800_50_16_1x  <-read.csv("res_800_50_16_1x_M4.csv")

df5_800_400_2_1x  <-read.csv("res_800_400_2_1x_M5.csv")
df5_800_200_4_1x  <-read.csv("res_800_200_4_1x_M5.csv")
df5_800_100_8_1x  <-read.csv("res_800_100_8_1x_M5.csv")
df5_800_50_16_1x  <-read.csv("res_800_50_16_1x_M5.csv")

df <- gdata::combine(
  df1_800_400_2_1x, 
  df1_800_200_4_1x, 
  df1_800_100_8_1x,
  df1_800_50_16_1x,
  
  df4_800_400_2_1x, 
  df4_800_200_4_1x, 
  df4_800_100_8_1x,
  df4_800_50_16_1x,
  
  df5_800_400_2_1x, 
  df5_800_200_4_1x, 
  df5_800_100_8_1x,
  df5_800_50_16_1x
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
df$RBias <- ifelse(df$Parameter=="psi", (df$V3-Mpsi)/Mpsi, df$RBias)
df$RBias <- ifelse(df$Parameter=="Sigma2_phi", ((df$V3-Vphi)/Vphi),df$RBias)
df$RBias <- ifelse(df$Parameter=="cov_int_phi", ((df$V3-Cov_alpha.phi)/Cov_alpha.phi),df$RBias)

df <- rename(df, Q=source)
df$Q <- gsub('df','',df$Q)
df$Q <- as.factor(df$Q)
df[c('Model','N','ind','n_per_ind','repeats')] <- str_split_fixed(df$Q, '_', 5)
df <- df %>%  mutate(ind = fct_relevel(ind, "400","200","100","50"))

# MAE graph
df1<-na.omit(df)
summary <- df1 %>% group_by(Model,N,ind,n_per_ind,repeats,Parameter) %>%
  summarise(MAE = mean(abs(RBias)))

summary1 <- filter(summary, Parameter == "psi")
summary2 <- filter(summary, Parameter == "Sigma2_phi" )
summary3 <- filter(summary, Parameter == "cov_int_phi")

summary1$ind <- as.numeric(as.factor(summary1$ind))

p1<-ggplot(summary1,(aes(x=ind, y=MAE*100, col=Parameter, group=Model))) +
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=MAE*100,linetype=Model), size=0.7) +
  scale_x_continuous(name="Number of individuals/\nsocial partners", breaks=(1:4),labels=c("400\n2","200\n4","100\n8","50\n16"))+
  scale_y_continuous(name="Mean absolute error (%)", breaks = c(0,10,20,30,40), limits= c(0,44)) +
  scale_colour_manual(values=c("#E5C494"), labels=c("\U1D713")) +
  
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),         
        axis.text = element_text(color="black"),
        axis.title.x=element_blank(),
        legend.position = "none",
        legend.text=element_text(size=14),
        legend.title = element_text())


summary2$ind <- as.numeric(as.factor(summary2$ind))

p2<-ggplot(summary2,(aes(x=ind, y=MAE*100, col=Parameter, group=Model))) +
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=MAE*100,linetype=Model), size=0.7) +
  scale_x_continuous(name="Number of individuals/\nsocial partners", breaks=(1:4),labels=c("400\n2","200\n4","100\n8","50\n16"))+
  scale_y_continuous(name="Mean absolute error (%)", breaks = c(0,10,20,30,40), limits= c(0,44)) +
  scale_colour_manual(values=c("#8DA0CB"), labels=c("\U1D449\U1D719"),
                      guide="none") +
  scale_shape_discrete(labels=c("novel I&R", "Hybrid+EIV", "Hybrid+RS"))+
  scale_linetype_discrete(labels=c("novel I&R", "Hybrid+EIV", "Hybrid+RS"))+
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


summary3$ind <- as.numeric(as.factor(summary3$ind))

p3<-ggplot(summary3,(aes(x=ind, y=MAE*100, col=Parameter, group=Model))) +
  geom_point(aes(shape=Model), size=2)+
  geom_line(aes(y=MAE*100,linetype=Model), size=0.7) +
  scale_x_continuous(name="Number of individuals/\nsocial partners", breaks=(1:4),labels=c("400\n2","200\n4","100\n8","50\n16"))+
  scale_y_continuous(name="Mean absolute error (%)", breaks = c(0,10,20,30,40), limits= c(0,44)) +
  scale_colour_manual(values=c("#E78AC3"), labels=c("\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D719\U2769")) +
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
        legend.text=element_text(size=14),
        legend.title = element_text())


## Boxplot ----
df1 <- filter(df, Parameter=="psi")
df1$ind <- as.numeric(as.factor(df1$ind))

plot1<-ggplot(df1,aes(x = ind, y = V3, fill=Model,col= Parameter, group=interaction(Model,ind,Parameter))) + 
  geom_boxplot(size=0.8, outlier.shape = NA) +
  scale_x_continuous(name="Number of individuals/\nsocial partners", breaks=c(1:4),labels=c("400\n2","200\n4","100\n8","50\n16"))+
  scale_y_continuous(name="Bias and precision", breaks = c(0,0.1,0.2,0.3,0.4), limits= c(0,0.44)) +
  geom_hline(aes(yintercept = Mpsi), linetype = 'dashed', col = "#E5C494") +
  scale_colour_manual(values="#E5C494", labels=c("\U1D713")) +
  scale_fill_brewer(palette="Greys") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),         
    axis.text = element_text(color="black"),
    axis.title.x=element_blank(),
    legend.position = "none",
    legend.text=element_text(size=14),
    legend.title = element_text())


df2 <- filter(df, Parameter == "Sigma2_phi")
df2$ind <- as.numeric(as.factor(df2$ind))

plot2<-ggplot(df2,aes(x = ind, y = V3, fill=Model,col= Parameter, group=interaction(Model,ind,Parameter))) + 
  geom_boxplot(size=0.8, outlier.shape = NA) +
  scale_x_continuous(name="Number of individuals/\nsocial partners", breaks=c(1:4),labels=c("400\n2","200\n4","100\n8","50\n16"))+
  scale_y_continuous(name="Bias and precision", breaks = c(0,0.1,0.2), limits= c(-0.02,0.25)) +
  geom_hline(aes(yintercept = Vphi), linetype = 'dashed', col = "#8da0cb") +
  scale_colour_manual(values=c("#8DA0CB"), labels=c("\U1D449\U1D719"),
                      guide="none") +
  scale_fill_brewer(palette = "Greys", guide="none") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),         
    axis.text = element_text(color="black"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position = "right",
    legend.text=element_text(size=14),
    legend.title = element_text())


df3 <- filter(df, Parameter == "cov_int_phi")
df3$ind <- as.numeric(as.factor(df3$ind))

plot3<-ggplot(df3,aes(x = ind, y = V3, col= Parameter,fill=Model, group=interaction(Model,ind,Parameter))) + 
  geom_boxplot(size=0.8, outlier.shape = NA) +
  scale_x_continuous(name="Number of individuals/\nsocial partners", breaks=c(1:4),labels=c("400\n2","200\n4","100\n8","50\n16"))+
  scale_y_continuous(name="Bias and precision", breaks = c(0,0.1,0.2), limits= c(-0.02,0.25)) +
  geom_hline(aes(yintercept = Cov_alpha.phi), linetype = 'dashed', col =  "#E78AC3") +
  labs(y = "Bias and precision", x = "Sample size") +
  scale_colour_manual(values=c("#E78AC3"), labels=c("\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D719\U2769"),
                      guide="none") +
  scale_fill_brewer(palette = "Greys", labels=c("novel I&R", "Hybrid+EIV", "Hybrid+RS")) +
  theme_bw(base_size = 10)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(color="black"),
    axis.text = element_text(color="black"),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position = "right",
    legend.text=element_text(size=10),
    legend.title = element_text())


patch1 <- (p1+p2+p3+plot1 + plot2+plot3) + plot_layout(guides = 'collect')
wrap_elements(panel = patch1) + 
  labs(tag = "Number of individuals/\nsocial partners") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  ) 


