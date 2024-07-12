library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(forcats)

setwd("~/Social-impact-and-responsivenesss/Sample_size")

df_400_100_4_1x <-read.csv("res_400_100_4_1x.csv")
df_800_200_4_1x <-read.csv("res_800_200_4_1x.csv")
df_1600_400_4_1x<-read.csv("res_1600_400_4_1x.csv")
df_3200_800_4_1x<-read.csv("res_3200_800_4_1x.csv")
df_6400_1600_4_1x<-read.csv("res_6400_1600_4_1x.csv")

df <- gdata::combine(
  df_400_100_4_1x,
  df_800_200_4_1x,
  df_1600_400_4_1x,
  df_3200_800_4_1x,
  df_6400_1600_4_1x
  )

Mad<- 1
Mpsi <- 0.1

df <- rename(df, Q=source)
df$Q <- gsub('df','',df$Q)
df$Q <- as.factor(df$Q)
df[c('Q','N','ind','n_per_ind','repeats')] <- str_split_fixed(df$Q, '_', 5)
df <- df %>%  mutate(N = fct_relevel(N, "400","800","1600","3200","6400"))
df$N <- as.numeric(as.factor(df$N))

df1 <- filter(df, Parameter=="B_0"|Parameter=="psi")
df2 <- filter(df, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi")
df3 <- filter(df, Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")

p1<-ggplot(df1,aes(x=N, y=V3, fill=Parameter, group=interaction(N,Parameter))) +
  geom_boxplot(size=0.6,outlier.shape=NA) +
  scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1,1.2)) +
  scale_x_continuous(name="Total number of observations", breaks=c(1:5), labels=c(400,800,1600,3200,6400)) +
  geom_hline(aes(yintercept = Mad), linetype='dashed', col="#B3B3B3") +
  geom_hline(aes(yintercept = Mpsi), linetype='dashed', col="#E5C494") +
  labs(y="Bias and precision", x="Sample size") +
  scale_fill_manual(values = c("#B3B3B3","#E5C494"), labels=c("\U1D6FD","\U1D713")) +
  theme_bw(base_size=10) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.title.x=element_blank(),
        axis.text=element_text(color="black"),
        legend.position = "right",
        legend.text=element_text(size=10),
        legend.title = element_text())

p2<-ggplot(df2,aes(x=N, y=V3, fill=Parameter, group=interaction(N,Parameter))) +
  geom_boxplot(size=0.6,outlier.shape=NA) +
  scale_x_continuous(name="Total number of observations", breaks=c(1:5), labels=c(400,800,1600,3200,6400))+
  scale_y_continuous(limits= c(0.0,0.33)) +
  geom_hline(aes(yintercept = Vad), linetype = 'dashed', col = "#66c2a5") +
  geom_hline(aes(yintercept = Vpsi+0.001), linetype = 'dashed', col = "#fc8d62") +
  geom_hline(aes(yintercept = Mpsi^2 * Vx + Vas-0.001), linetype = 'dashed', col = "#8da0cb") +
  scale_fill_manual(values=c("#66C2A5", "#8DA0CB", "#FC8D62"), labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
  theme_bw(base_size=10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.y=element_blank(),
    axis.line = element_line(),
    axis.title.x=element_blank(),
    axis.text=element_text(color="black"),
    axis.ticks = element_line(),
    legend.position = "right",
    legend.text=element_text(size=14),
    legend.title = element_text())

p3<-ggplot(df3,aes(x=N, y=V3, fill=Parameter, group=interaction(N,Parameter))) +
  geom_boxplot(size=0.6,outlier.shape=NA) +
  scale_x_continuous(name="Total number of observations", breaks=c(1:5), labels=c(400,800,1600,3200,6400))+
  scale_y_continuous(limits= c(-0.15,0.15))+
  geom_hline(aes(yintercept = Cdpsi), linetype = 'dashed', col ="#A6D854") +
  geom_hline(aes(yintercept = Cds+Mpsi*Cdx), linetype = 'dashed', col =  "#E78AC3") +
  geom_hline(aes(yintercept = Cspsi+Mpsi*Cpsix), linetype = 'dashed', col = "#FFD92F") +
  labs(y = "Bias and precision", x = "Sample size") +
  scale_fill_manual(values=c("#E78AC3","#A6D854", "#FFD92F"), labels=
                        c("\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D719\U2769",
                          "\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D713\U2769",
                          "\U1D436\U1D45C\U1D463\U2768\U1D713\U002C\U1D719\U2769")) +
  theme_bw(base_size=10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    axis.title.x=element_blank(),
    axis.text=element_text(color="black"),
    axis.title.y=element_blank(),
    legend.position = "right",
    legend.text=element_text(size=14),
    legend.title = element_text())

patch1 <- (p1 + p2 + p3) + plot_layout(guides = 'collect')
wrap_elements(panel = patch1) +
  labs(tag = "Total number of observations") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  ) 
