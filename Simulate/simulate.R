
require(MASS)

# Means
Malpha   = 1      # mean focal effect
Mepsilon = 0      # mean partner effect
Mpsi     = 0.3    # interaction coefficient psi
Mx       = 0      # mean-centered covariate/opponent trait

# Variances
Valpha   = 0.2    # variance direct effect (mean behaviour)
Vepsilon = 0.01   # variance social partner effect (residual impact) (so that total impact sums to 0.1)
Vpsi     = 0.1    # variance in slopes (social responsiveness)
Vx       = 1      # variance social partner phenotype (impact covariate)
Ve       = 0.6    # residual variance
Vxe      = 0.1    # measurement error or variation of labile trait

# Correlations between variance components
r_alpha_epsilon  =  0       # part of Cov int-phi
r_alpha_psi      =  -0.6    # part of Cov int-psi 
r_epsilon_psi    =  -0.6    # part of Cov psi-phi 
r_alpha_x        =  0.6     # part of Cov int-phi 
r_psi_x          =  -0.6    # part of Cov psi-phi 
r_epsilon_x      =  0



# Random balanced design function
sampling_design_balanced <- function(n_ind,partners,repeats){
  ind_seq <- 1:n_ind 
  IDi <- rep(ind_seq, partners) 
  IDj = vector("list", length = n_ind*partners) 
  for (i in 1:partners) {
    seq_j <- c(ind_seq[-c(1:i)], ind_seq[1:i])
    IDj[[i]] <- seq_j   }
  IDj = do.call(c, IDj)
  df <- data.frame(IDi,IDj)
  df <- df[rep(seq_len(nrow(df)), repeats),] 
  df
}

# Means & Variance G-matrix for multivariate normal distribution
M<-c(Malpha,Mepsilon,Mpsi,Mx)
G<-matrix(NA,4,4)

G[1,1]<-Valpha
G[2,2]<-Vepsilon
G[3,3]<-Vpsi
G[4,4]<-Vx

G[1,2]<-G[2,1]<-r_alpha_epsilon*(sqrt(Valpha*Vepsilon)) 
G[1,3]<-G[3,1]<-r_alpha_psi*(sqrt(Valpha*Vpsi)) 
G[2,3]<-G[3,2]<-r_epsilon_psi*(sqrt(Vepsilon*Vpsi))
G[4,1]<-G[1,4]<-r_alpha_x*(sqrt(Valpha*Vx))
G[4,2]<-G[2,4]<-r_epsilon_x*(sqrt(Vepsilon*Vx))
G[4,3]<-G[3,4]<-r_psi_x*(sqrt(Vpsi*Vx))

# Create phenotypic traits of individuals and outcome of interactions (requires study design)
sim<-function(df, M, G, Ve) {
  a=as.data.frame(mvrnorm(n_ind, M, G))
  colnames(a)<-c("alpha", "epsilon", "psi", "x")
  a$ID=1:n_ind
  df$alpha_i=a[match(df$IDi,a$ID),"alpha"]
  df$epsilon_j=a[match(df$IDj,a$ID),"epsilon"]
  df$psi_i=a[match(df$IDi,a$ID),"psi"]
  df$xj=a[match(df$IDj,a$ID),"x"]
  df$x_ijk <- df$xj + rnorm(nrow(df),0, sqrt(Vxe))
  df$e_ijk<-rnorm(nrow(df),0, sqrt(Ve))
  
  df$z_i=df$alpha_i + df$psi_i*df$xj + df$epsilon_j + df$e_ijk # Phenotypic Equation
  df
}


n_ind     = 200 # number of individuals 
partners  = 2   # number of partners per individual
repeats   = 8   # number of repeats of unique dyadic interaction
iterations= 1000   # number of datasets


# Simulate datasets and create lists
df<-sampling_design_balanced(n_ind, partners, repeats) # see Random balanced design function
dfl <- list()
for(i in 1:iterations){
  df2<-sim(df, M, G, Ve)
  xi<-aggregate(df2$x_ijk, list(IDi = df2$IDj), mean)
  dfl[[i]]<-list(n_obs = nrow(df2),
                 n_ind=length(unique((df2$IDi))),
                 individual=df2$IDi,
                 opponent=df2$IDj,
                 xj=df2$x_ijk,  # model input is the observed xj
                 xi=xi$x,
                 z=df2$z_i)
}

# List saving: dfl_N(total observations)_N(indidviduals)_N(social partners)_N(repeats)
saveRDS(
  dfl,
  file = paste0(
    "Data/dfl_",
    n_ind*partners*repeats,"_",n_ind, "_", partners, "_", repeats, "x.RDS"
  )
)
