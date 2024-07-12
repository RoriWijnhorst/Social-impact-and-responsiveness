library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(forcats)

df1_800_400_2_1x  <-read.csv("res_800_400_2_1x.csv")
df1_800_200_4_1x  <-read.csv("res_800_200_4_1x.csv")
df1_800_100_8_1x  <-read.csv("res_800_100_8_1x.csv")
df1_800_50_16_1x  <-read.csv("res_800_50_16_1x.csv")

df2_800_200_4_1x <-read.csv("res_800_200_4_1x.csv")
df2_800_100_4_2x <-read.csv("res_800_100_4_2x.csv")
df2_800_50_4_4x  <-read.csv("res_800_50_4_4x.csv")

df3_800_100_8_1x <-read.csv("res_800_100_8_1x.csv")
df3_800_100_4_2x <-read.csv("res_800_100_4_2x.csv")
df3_800_100_2_4x <-read.csv("res_800_100_2_4x.csv")

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

df$Q <- df$source
df$Q <- gsub('df','',df$Q)
df$Q <- as.factor(df$Q)
df[c('Q','N','ind','n_per_ind','repeats')] <- str_split_fixed(df$Q, '_', 5)
df <- df %>%  mutate(ind = fct_relevel(ind, "400","200","100","50"))

## boxplots ----
dfv_Q1 <- filter(df, Q==1, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi")
dfv_Q1$ind <- as.numeric(as.factor(dfv_Q1$ind))
MAE <- dfv_Q1 %>% group_by(ind,Parameter) %>% summarize(MAE = mean(abs(RBias))) %>%   ungroup() 


plot1<-ggplot(dfv_Q1, aes(x = ind, y = V3, fill = Parameter, group=interaction(ind,Parameter))) + 
  geom_boxplot(size=0.5, outlier.shape = NA) +
  scale_x_continuous(name="Number of individuals/\nsocial partners", breaks=c(1:4), labels=c("400\n2","200\n4","100\n8","50\n16"))+
  scale_y_continuous(name="Bias and precision", breaks = c(0,0.1,0.2,0.3), limits= c(0,0.35)) +
  geom_hline(aes(yintercept = Vad), linetype = 'longdash', col = "#66c2a5") +
  geom_hline(aes(yintercept = Vpsi+0.001), linetype = 'longdash', col = "#fc8d62") +
  geom_hline(aes(yintercept = Mpsi^2 * Vx + Vas-0.001), linetype = 'longdash', col = "#8da0cb") +
  scale_fill_manual(values=c("#66C2A5", "#8DA0CB", "#FC8D62"), labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
  annotate("text", x=2.5, y=0.35, label= "repeats = 1x",size=4) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x=element_blank(),
    axis.text = element_text(color="black"),
    axis.text.x = element_blank(),
    axis.ticks = element_line(),
    legend.position = "none",
    legend.text=element_text(size=14),
    legend.title = element_text())

dfv_Q2 <- filter(df, Q==2, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi")
dfv_Q2$ind <- as.numeric(as.factor(dfv_Q2$ind))

plot2<-ggplot(dfv_Q2, aes(x = ind, y = V3, fill = Parameter, group=interaction(ind,Parameter))) + 
  geom_boxplot(size=0.5, outlier.shape = NA) +
  scale_x_continuous(name="Number of individuals/\nrepeats", breaks=c(2,3,4), labels=c("200\n1x","100\n2x","50\n4x"))+
  scale_y_continuous(name="Bias and precision", breaks = c(0,0.1,0.2,0.3), limits= c(0,0.35)) +
  geom_hline(aes(yintercept = Vad), linetype = 'longdash', col = "#66c2a5") +
  geom_hline(aes(yintercept = Vpsi+0.001), linetype = 'longdash', col = "#fc8d62") +
  geom_hline(aes(yintercept = Mpsi^2 * Vx + Vas-0.001), linetype = 'longdash', col = "#8da0cb") +
  scale_fill_manual(values=c("#66C2A5", "#8DA0CB", "#FC8D62"), labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
  annotate("text", x=3, y=0.35, label= "partners = 4",size=4) +
  labs(title=expression(''~italic(N)~'= 800'),)+
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text = element_text(color="black"),
    axis.text.x = element_blank(),
    axis.ticks = element_line(),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    legend.text=element_text(size=14),
    legend.title = element_text())

dfv_Q3 <- filter(df, Q==3, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi")
dfv_Q3$n_per_ind <- fct_relevel(dfv_Q3$n_per_ind, "8","4","2")
dfv_Q3$n_per_ind <- as.numeric(as.factor(dfv_Q3$n_per_ind))

plot3<-ggplot(dfv_Q3, aes(x = n_per_ind, y = V3, fill = Parameter, group=interaction(n_per_ind,Parameter))) + 
  geom_boxplot(size=0.5, outlier.shape = NA) +
  scale_x_continuous(name="Number of social partners/\nrepeats", breaks=c(1,2,3), labels=c("8\n1x","4\n2x","2\n4x"))+
  scale_y_continuous(name="Bias and precision", breaks = c(0,0.1,0.2,0.3), limits= c(0,0.35)) +
  geom_hline(aes(yintercept = Vad), linetype = 'longdash', col = "#66c2a5") +
  geom_hline(aes(yintercept = Vpsi+0.001), linetype = 'longdash', col = "#fc8d62") +
  geom_hline(aes(yintercept = Mpsi^2 * Vx + Vas-0.001), linetype = 'longdash', col = "#8da0cb") +
  scale_fill_manual(values=c("#66C2A5", "#8DA0CB", "#FC8D62"), labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
  annotate("text", x=2, y=0.35, label= "individuals = 100",size=4) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(color="black"),
    axis.text.x = element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.ticks = element_line(),
    legend.position = "right",
    legend.text=element_text(size=14),
    legend.title = element_text())


dfc_Q1 <- filter(df, Q==1, Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")
dfc_Q1$ind <- as.numeric(as.factor(dfc_Q1$ind))

p1<-ggplot(dfc_Q1,aes(x = ind, y = V3, fill = Parameter, group=interaction(ind,Parameter))) + 
  geom_boxplot(size=0.5, outlier.shape = NA) +
  scale_x_continuous(name="Number of individuals/\nsocial partners", breaks=c(1:4),labels=c("400\n2","200\n4","100\n8","50\n16"))+
  scale_y_continuous(name="Bias and precision", breaks = c(-0.1,0,0.1), limits= c(-0.14,0.14)) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed',  col="grey") +
  geom_hline(aes(yintercept = Cdpsi), linetype = 'longdash', col ="#A6D854") +
  geom_hline(aes(yintercept = Cds+Mpsi*Cdx), linetype = 'longdash', col =  "#E78AC3") +
  geom_hline(aes(yintercept = Cspsi+Mpsi*Cpsix), linetype = 'longdash', col = "#FFD92F") +
  labs(y = "Bias and precision", x = "Sample size") +
  scale_fill_manual(values=c("#E78AC3","#A6D854", "#FFD92F"), labels=
                      c("\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D719\U2769",
                        "\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D713\U2769",
                        "\U1D436\U1D45C\U1D463\U2768\U1D713\U002C\U1D719\U2769")) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.ticks = element_line(),
    axis.text = element_text(color="black"),
    legend.position = "none",
    legend.text=element_text(size=14),
    legend.title = element_text())

dfc_Q2 <- filter(df, Q==2, Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")
dfc_Q2$ind <- as.numeric(as.factor(dfc_Q2$ind))

p2<-ggplot(dfc_Q2,aes(x = ind, y = V3, fill = Parameter, group=interaction(ind,Parameter))) + 
  geom_boxplot(size=0.5, outlier.shape = NA) +
  scale_x_continuous(name="Number of individuals/\nrepeats", breaks=c(2,3,4), labels=c("200\n1x","100\n2x","50\n4x"))+
  scale_y_continuous(name="Bias and precision", breaks = c(-0.1,0,0.1), limits= c(-0.14,0.14)) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed',  col="lightgrey") +
  geom_hline(aes(yintercept = Cdpsi), linetype = 'longdash', col ="#A6D854") +
  geom_hline(aes(yintercept = Cds+Mpsi*Cdx), linetype = 'longdash', col =  "#E78AC3") +
  geom_hline(aes(yintercept = Cspsi+Mpsi*Cpsix), linetype = 'longdash', col = "#FFD92F") +
  labs(y = "Bias and precision", x = "Sample size") +
  scale_fill_manual(values=c("#E78AC3","#A6D854", "#FFD92F"), labels=
                      c("\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D719\U2769",
                        "\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D713\U2769",
                        "\U1D436\U1D45C\U1D463\U2768\U1D713\U002C\U1D719\U2769")) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    axis.line = element_line(),
    axis.title.y=element_blank(),
    axis.ticks = element_line(),
    axis.text = element_text(color="black"),
    legend.position = "none",
    legend.text=element_text(size=14),
    legend.title = element_text())


dfc_Q3 <- filter(df, Q==3, Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")
dfc_Q3$n_per_ind <- fct_relevel(dfc_Q3$n_per_ind, "8","4","2")
dfc_Q3$n_per_ind <- as.numeric(as.factor(dfc_Q3$n_per_ind))

p3<-ggplot(dfc_Q3,aes(x = n_per_ind, y = V3, fill = Parameter, group=interaction(n_per_ind,Parameter))) + 
  geom_boxplot(size=0.5, outlier.shape = NA) +
  scale_x_continuous(name="Number of social partners/\nrepeats", breaks=c(1,2,3), labels=c("8\n1x","4\n2x","2\n4x")) +
  scale_y_continuous(name="Bias and precision", breaks = c(-0.1,0,0.1), limits= c(-0.14,0.14)) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed',  col="lightgrey") +
  geom_hline(aes(yintercept = Cdpsi), linetype = 'longdash', col ="#A6D854") +
  geom_hline(aes(yintercept = Cds+Mpsi*Cdx), linetype = 'longdash', col =  "#E78AC3") +
  geom_hline(aes(yintercept = Cspsi+Mpsi*Cpsix), linetype = 'longdash', col = "#FFD92F") +
  labs(y = "Bias and precision", x = "Sample size") +
  scale_fill_manual(values=c("#E78AC3","#A6D854", "#FFD92F"), labels=
                      c("\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D719\U2769",
                        "\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D713\U2769",
                        "\U1D436\U1D45C\U1D463\U2768\U1D713\U002C\U1D719\U2769")) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.title.y=element_blank(),
    axis.text = element_text(color="black"),
    axis.ticks = element_line(),
    legend.position = "right",
    legend.text=element_text(size=14),
    legend.title = element_text())


patch1 <- (plot1+plot2+plot3)/(p1 + p2 + p3) + plot_layout(guides = 'collect')
wrap_elements(panel = patch1) + 
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  ) 




