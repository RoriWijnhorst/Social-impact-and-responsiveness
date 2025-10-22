library(stringr)
library(patchwork)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(forcats)
getwd()
df_400_50_8_1x <-read.csv("Output/res_400_50_8_1x_M1.csv")
df_800_100_8_1x <-read.csv("Output/res_800_100_8_1x_M1.csv")
df_1600_200_8_1x<-read.csv("Output/res_1600_200_8_1x_M1.csv")
df_3200_400_8_1x<-read.csv("Output/res_3200_400_8_1x_M1.csv")
df_6400_800_8_1x<-read.csv("Output/res_6400_800_8_1x_M1.csv")
df <- gdata::combine(
    df_400_50_8_1x,
    df_800_100_8_1x,
    df_1600_200_8_1x,
    df_3200_400_8_1x,
     df_6400_800_8_1x
   )

# Means
 Malpha   = 1      # mean focal effect
 Mepsilon = 0      # mean partner effect
 Mpsi     = 0.3   # interaction coefficient psi
 Mx       = 0      # mean-centered covariate/opponent trait
 # Variances
 Valpha   = 0.2    # variance direct effect (mean behaviour)
 Vepsilon = 0.01  # variance social partner effect (residual impact)
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
       df[c('Q','N','ind','n_per_ind','repeats')] <- str_split_fixed(df$Q, '_', 5)
       df <- df %>%  mutate(ind = fct_relevel(ind,"50", "100","200","400","800"))

       
     
       Vphi <- Mpsi^2 * Vx + Vepsilon
       Cov_alpha.phi <- r_alpha_epsilon*(sqrt(Valpha*Vepsilon))+ Mpsi*r_alpha_x*(sqrt(Valpha*Vx))
       Vphi <- Mpsi^2 * Vx + Vepsilon
       Cov_alpha.phi <- r_alpha_epsilon*(sqrt(Valpha*Vepsilon))+ Mpsi*r_alpha_x*(sqrt(Valpha*Vx))
       Cov_psi.phi <- r_epsilon_psi*(sqrt(Vepsilon*Vpsi)) + Mpsi*r_psi_x*(sqrt(Vpsi*Vx))
       Cov_alpha.psi <- r_alpha_psi*(sqrt(Valpha*Vpsi))
       df$RBias <- NA
       df$RBias <- ifelse(df$Parameter=="B_0", (df$X50.-Malpha)/Malpha, df$RBias)
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
      
       RMAD <- function(x) {
         mean_x <- mean(x)
         mae <- mean(abs(x - mean_x))
         mean_abs <- mean(abs(x))
         
         if (mean_abs == 0) return(NA_real_)
         return(mae / mean_abs)
       }
       # Compute summary statistics
       summary <- df %>% group_by(Q,N,ind,n_per_ind,repeats, Parameter) %>%
           summarise(
             Rhat=median(Rhat),
             N_eff=median(n_eff),
             mean_bias = mean(RBias), # Single point
          sd_bias = sd(RBias),# Dispersion measure
               RMAD = RMAD(X50.),
               .groups = 'drop'
             )
       
       
       

       
       # Plot bias and dispersion
       df_Q1 <- filter(summary, Parameter=="psi"|Parameter=="B_0")
       df_Q1$ind <- as.numeric(as.factor(df_Q1$ind))
   
      # Plot bias and dispersion
         dfv_Q1 <- filter(summary, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi")
      dfv_Q1$ind <- as.numeric(as.factor(dfv_Q1$ind))
       # Plot bias and dispersion
       dfc_Q1 <- filter(summary, Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")
       dfc_Q1$ind <- as.numeric(as.factor(dfc_Q1$ind))
       
       
       # === Bias Plot ===
       a <- ggplot(df_Q1, aes(x = as.factor(ind), y = mean_bias*100, color = Parameter, group = Parameter)) + 
         geom_hline(yintercept = 0, linetype = "dashed", color = "black") +    
         geom_line(size=0.7) +         
         geom_point(size=2) +        +
         scale_x_discrete(name = "Total observations", labels=c("400","800","1600","3200","6400")) +
         scale_y_continuous(name = "Bias (%)", limits = c(-5, 10), breaks = seq(-5, 10, 5)) +
         scale_color_manual(values = c("#B3B3B3","#E5C494"), labels=c("\U1D6FD","\U1D713")) +
         theme_bw(base_size=10) +
         theme(
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           axis.title.x = element_blank(),
           axis.line.x = element_line(color="black"), 
           axis.line.y = element_line(color="black"),  
           axis.text = element_text(color="black"),
           axis.text.x = element_blank(),
           axis.ticks = element_line(color="black"),
           legend.position = "right",
           legend.text = element_text(size=14),
           legend.title = element_text()
         )
       
       
       b <- ggplot(dfv_Q1, aes(x = as.factor(ind), y = mean_bias*100, color = Parameter, group = Parameter)) + 
           geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
           geom_line(size=0.7) +         
           geom_point(size=2) +        +
           scale_x_discrete(name = "Total observations", labels=c("400","800","1600","3200","6400")) +
           scale_y_continuous(name = "Bias (%)", limits = c(-5, 10), breaks = seq(-5, 10, 5)) +
         scale_color_manual(values=c("#2E6F5D", "#66C2A5", "#B2E1D3"), 
                            labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
           theme_bw(base_size = 10) +
         theme(
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           axis.title.y = element_blank(),
           axis.line.x = element_line(color="black"), 
           axis.line.y = element_line(color="black"),  
           axis.title = element_blank(),  
           axis.text = element_text(color="black"),
           axis.text.x = element_blank(),
           axis.ticks = element_line(color="black"),
           legend.position = "right",
           legend.text = element_text(size=14),
           legend.title = element_text()
         )
      
       
      
       c<-ggplot(dfc_Q1, aes(x = as.factor(ind), y = mean_bias*100, color = Parameter, group = Parameter)) + 
               geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
           geom_line(size=0.7) +         
            geom_point(size=2) +        +
         scale_x_discrete(name = "Total observations", labels=c("400","800","1600","3200","6400")) +
         scale_y_continuous(name="Bias (%)", limits=c(-35,5), breaks = seq(-35,5,5))+
         scale_color_manual(values=c("#3F4176","#6D6EB0", "#C0C0DC"), labels=
                             c("\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D719\U2769",
                               "\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D713\U2769",
                               "\U1D436\U1D45C\U1D463\U2768\U1D713\U002C\U1D719\U2769")) +
         theme_bw(base_size = 10) +
         theme(
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           axis.title.y = element_blank(),
           axis.line.x = element_line(color="black"), 
           axis.line.y = element_line(color="black"),  
           axis.title = element_blank(),  
           axis.text = element_text(color="black"),
           axis.text.x = element_blank(),
           axis.ticks = element_line(color="black"),
           legend.position = "right",
           legend.text = element_text(size=14),
           legend.title = element_text()
         )
       
       
       df1 <- filter(df, Parameter=="B_0"|Parameter=="psi")
       df2 <- filter(df, Parameter=="Sigma2_intercept"|Parameter=="Sigma2_psi"|Parameter=="Sigma2_phi")
       df3 <- filter(df, Parameter=="cov_int_psi"|Parameter=="cov_int_phi"|Parameter=="cov_psi_phi")
       
       p1<-ggplot(df1,aes(x=as.numeric(ind), y=X50., fill=Parameter, group=interaction(N,Parameter))) +
         geom_hline(aes(yintercept = Malpha), linetype='dashed', col="#B3B3B3") +
         geom_hline(aes(yintercept = Mpsi), linetype='dashed', col="#E5C494") +
         geom_boxplot(size=0.4,outlier.shape=NA) +
         scale_y_continuous(limits=c(0,1.3), breaks= c(0,0.2,0.4,0.6,0.8,1,1.2)) +
         scale_x_continuous(name="Total number of observations", breaks=c(1:5), labels=c("50\n400","100\n800","200\n1600","400\n3200","800\n6400")) +
         labs(y="Model estimates", x="Sample size") +
         scale_fill_manual(values = c("#B3B3B3","#E5C494"), labels=c("\U1D6FD","\U1D713")) +
         theme_bw(base_size=10) +
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               axis.line = element_line(),
               axis.ticks = element_line(color="black"),
               axis.title.x=element_blank(),
               axis.text=element_text(color="black"),
               legend.position = "none",
               legend.text=element_text(size=10),
               legend.title = element_text())
       
       p2<-ggplot(df2,aes(x=as.numeric(ind), y=X50., fill=Parameter, group=interaction(N,Parameter))) +
         geom_hline(aes(yintercept = Valpha), linetype = 'dashed', col = "#2E6F5D") +
         geom_hline(aes(yintercept = Vpsi+0.001), linetype = 'dashed', col = "#66C2A5") +
         geom_hline(aes(yintercept = (Mpsi^2 * Vx + Vepsilon)-0.001), linetype = 'dashed', col = "#B2E1D3") +
         geom_boxplot(size=0.4,outlier.shape=NA) +
         scale_x_continuous(name="Total number of observations", breaks=c(1:5), labels=c("50\n400","100\n800","200\n1600","400\n3200","800\n6400"))+
         scale_y_continuous(limits= c(0.0,0.38)) +
        scale_fill_manual(values=c("#2E6F5D", "#66C2A5", "#B2E1D3"), labels=c("\U1D449\U1D6FC", "\U1D449\U1D719", "\U1D449\U1D713")) +
         theme_bw(base_size=10) +
         theme(
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           axis.title.y = element_blank(),
           axis.line.x = element_line(color="black"),   # ensure x-axis baseline
           axis.line.y = element_line(color="black"),   # ensure y-axis baseline
           axis.title.x = element_blank(),       # or element_blank() if you want no title
           axis.text = element_text(color="black"),
           axis.ticks = element_line(color="black"),
           legend.position = "none",
           legend.text = element_text(size=14),
           legend.title = element_text()
         )
       
       
       p3<-ggplot(df3,aes(x=as.numeric(ind), y=X50., fill=Parameter, group=interaction(N,Parameter))) +
         geom_hline(aes(yintercept = r_alpha_epsilon*(sqrt(Valpha*Vepsilon))+ Mpsi*r_alpha_x*(sqrt(Valpha*Vx))), 
                    linetype = 'dashed', col =  "#3F4176") +
         geom_hline(aes(yintercept = r_alpha_psi*(sqrt(Valpha*Vpsi))), 
                    linetype = 'dashed', col ="#6D6EB0") +
         geom_hline(aes(yintercept = r_epsilon_psi*(sqrt(Vepsilon*Vpsi))    + Mpsi*r_psi_x*(sqrt(Vpsi*Vx))), 
                    linetype = 'dashed', col =  "#C0C0DC") +
         geom_boxplot(size=0.4,outlier.shape=NA) +
         scale_x_continuous(name="Total number of observations", breaks=c(1:5), labels=c("50\n400","100\n800","200\n1600","400\n3200","800\n6400"))+
         scale_y_continuous(limits= c(-0.18,0.18))+
          labs(y = "Bias and precision", x = "Sample size") +
         scale_fill_manual(values=c("#3F4176","#6D6EB0", "#C0C0DC"), labels=
                             c("\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D719\U2769",
                               "\U1D436\U1D45C\U1D463\U2768\U1D6FC\U002C\U1D713\U2769",
                               "\U1D436\U1D45C\U1D463\U2768\U1D713\U002C\U1D719\U2769")) +
         theme_bw(base_size=10) +
         theme(
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           axis.line = element_line(),
           axis.ticks = element_line(color="black"),
           axis.title.x=element_blank(),
           axis.text=element_text(color="black"),
           axis.title.y=element_blank(),
           legend.position = "none",
           legend.text=element_text(size=14),
           legend.title = element_text())
       p3
       
       patch1 <- (a + b + c) / (p1 + p2 + p3) + plot_layout(guides = 'collect')
       wrap_elements(panel = patch1) +
         labs(tag = "Number of individuals/\nTotal observations") +
         theme(
           plot.tag = element_text(size = rel(1)),
           plot.tag.position = "bottom"
         )
       
       