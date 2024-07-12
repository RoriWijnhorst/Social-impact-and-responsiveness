library(parallel)
library(dplyr)
numCores <- 6
cl <- makeCluster(numCores)
clusterEvalQ(cl, {
  library(rstan)
  library(dplyr)
})

# M1 ----------------------
write( "data {
//Phenotype dataset

  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or aantal rows
   int<lower=0> n_ind; //number of individuals
   
  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent[n_obs];  //  Individual ID opponent repeated obs
   
   
  // Continuous variables
    real xj[n_obs];  // observed covariate (opponent trait) measured with error
     real z[n_obs];  // phenotypic observations
 }
 
 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
     real B_x; //intercept covariate model
     real B_0; //intercept
     real psi; //slope
     
     
   // Random effects
   matrix[4,n_ind]         zI; //(intercepts and slopes, res_impact for each individual)
   vector<lower=0>[4]      sigma_I; // sd  intercepts, slopes, res_impact
   cholesky_factor_corr[4] L;  // factor to estimate covariance int-slopes
   real<lower=0> sigma_e;
   real<lower=0> sigma_ex;
 }

 transformed parameters{
    matrix[4,n_ind] I; //  Unscaled blups intercept and slope and res_impact for each individual
    real e_z[n_obs]; // predicted values for phenotype
    real e_x[n_obs]; // predicted values for covariate
    I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value
    
     for (i in 1:n_obs) {
     e_x[i] = B_x +  I[4, focal[i]];
     }

   for (i in 1:n_obs) {
     e_z[i]  = B_0  +  I[1, individual[i]] + (psi + I[2, individual[i]])*I[4, opponent[i]] + I[3, opponent[i]];
   }
 }
 
model {
// Create vector of predicted values
  to_vector([B_0]) ~ normal(0, 1);
  to_vector([B_x]) ~ normal(0,1);
  to_vector([psi]) ~ normal(0, 1);
  

 // Random effects distribution
    to_vector(zI) ~ normal(0, 1);
    to_vector(sigma_I) ~ normal(0, 1);
    L ~ lkj_corr_cholesky(4);
    
    
 // Likelihood function
    for (i in 1:n_obs)
     z[i]~normal(e_z[i], sigma_e);
     
     for (i in 1:n_obs)
     xj[i]~normal(e_x[i], sigma_ex); 
     
 }

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_psi;
real<lower=0> Sigma2_epsilon;
real<lower=0> Sigma2_x;
real<lower=0> Sigma2_phi;
real<lower=0> Sigma2_e;
real<lower=0> Sigma2_ex;

real cov_1; //   psi-int
real cov_2; //   int-resimpact
real cov_3; //   psi-resimpact
real cov_4; //   x-int
real cov_5; //   x-psi
real cov_6; //   x-resimpact
real cov_int_psi;
real cov_int_phi;
real cov_psi_phi;

matrix[4, 4]  Omega_I;

Sigma2_intercept=sigma_I[1]^2; 
Sigma2_psi=sigma_I[2]^2;
Sigma2_epsilon=sigma_I[3]^2;
Sigma2_x=sigma_I[4]^2;

Sigma2_e=sigma_e^2;
Sigma2_ex=sigma_ex^2;

Omega_I = L * L';             
cov_1 = Omega_I[1,2]*sqrt(Sigma2_psi*Sigma2_intercept);   
cov_2 = Omega_I[1,3]*sqrt(Sigma2_epsilon*Sigma2_intercept);
cov_3 = Omega_I[2,3]*sqrt(Sigma2_psi*Sigma2_epsilon);
cov_4 = Omega_I[1,4]*sqrt(Sigma2_x*Sigma2_intercept);
cov_5 = Omega_I[2,4]*sqrt(Sigma2_x*Sigma2_psi);
cov_6 = Omega_I[3,4]*sqrt(Sigma2_x*Sigma2_epsilon);
cov_int_psi = cov_1;
cov_int_phi = cov_2 + psi*cov_4;  
cov_psi_phi = cov_3 + psi*cov_5; 
Sigma2_phi=sigma_I[3]^2 + psi^2*sigma_I[4]^2;  // variance total impact, also depends on covariance between x and epsilon

 }", file="RR.stan")


# data 300 ----------------------------
stan_func<-function(dfl_300_75_4_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_300_75_4_1x,
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept",  "Sigma2_psi", "Sigma2_x",
                      "Sigma2_epsilon", "Sigma2_phi", "cov_int_psi", "cov_int_phi",
                      "cov_psi_phi", "cov_1",  "cov_2", "cov_3", "cov_4", "cov_5",
                      "cov_6", "Sigma2_e", "Sigma2_ex"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_300_75_4_1x <- do.call(rbind,parSapply(cl, dfl_300_75_4_1x, stan_func))
res_300_75_4_1x <- as.data.frame(res_300_75_4_1x)
res_300_75_4_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept",  "Sigma2_psi", "Sigma2_x",
                              "Sigma2_epsilon", "Sigma2_phi", "cov_int_psi", "cov_int_phi",
                              "cov_psi_phi", "cov_1",  "cov_2", "cov_3", "cov_4", "cov_5",
                              "cov_6", "Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_300_75_4_1x, "res_300_75_4_1x_M1.csv")


## data 600 ----------------------------
stan_func<-function(dfl_600_150_4_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_600_150_4_1x,
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept",  "Sigma2_psi", "Sigma2_x",
                      "Sigma2_epsilon", "Sigma2_phi", "cov_int_psi", "cov_int_phi",
                      "cov_psi_phi", "cov_1",  "cov_2", "cov_3", "cov_4", "cov_5",
                      "cov_6", "Sigma2_e", "Sigma2_ex"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_600_150_4_1x <- do.call(rbind,parSapply(cl, dfl_600_150_4_1x, stan_func))
res_600_150_4_1x <- as.data.frame(res_600_150_4_1x)
res_600_150_4_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept",  "Sigma2_psi", "Sigma2_x",
                             "Sigma2_epsilon", "Sigma2_phi", "cov_int_psi", "cov_int_phi",
                             "cov_psi_phi", "cov_1",  "cov_2", "cov_3", "cov_4", "cov_5",
                             "cov_6", "Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_600_150_4_1x, "res_600_150_4_1x_M1.csv")


## data 1200 ----------------------------
stan_func<-function(dfl_1200_300_4_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_1200_300_4_1x, 
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept",  "Sigma2_psi", "Sigma2_x", 
                      "Sigma2_epsilon", "Sigma2_phi", "cov_int_psi", "cov_int_phi", 
                      "cov_psi_phi", "cov_1",  "cov_2", "cov_3", "cov_4", "cov_5", 
                      "cov_6", "Sigma2_e", "Sigma2_ex"), 
             control=list(max_treedepth = 15), chains = 1, iter = 5000,  
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_1200_300_4_1x <- do.call(rbind,parSapply(cl, dfl_1200_300_4_1x, stan_func)) 
res_1200_300_4_1x <- as.data.frame(res_1200_300_4_1x)
res_1200_300_4_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept",  "Sigma2_psi", "Sigma2_x", 
                             "Sigma2_epsilon", "Sigma2_phi", "cov_int_psi", "cov_int_phi", 
                             "cov_psi_phi", "cov_1",  "cov_2", "cov_3", "cov_4", "cov_5", 
                             "cov_6", "Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_1200_300_4_1x, "res_1200_300_4_1x_M1.csv")



## M2------------------------
write("data {
//Phenotype dataset

  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or aantal rows
   int<lower=0> n_ind; //number of individuals

  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent[n_obs];  //  Individual ID opponent repeated obs

  // Continuous variables
    real xj[n_obs];  // observed covariate (opponent trait) measured with error
     real z[n_obs];  // phenotypic observations
 }

 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
     real B_0; //intercept


   // Random effects
   matrix[2,n_ind]         zI; //(intercepts and slopes, res_impact for each individual)
   vector<lower=0>[2]      sigma_I; // sd  intercepts, slopes, res_impact
   cholesky_factor_corr[2] L;  // factor to estimate covariance int-slopes
   real<lower=0> sigma_e;
   real<lower=0> sigma_ex;
 }

 transformed parameters{
    matrix[2,n_ind] I; //  Unscaled blups intercept and slope and res_impact for each individual
    real e_z[n_obs]; // predicted values for covariate
    I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value


   for (i in 1:n_obs) {
     e_z[i]  = B_0  +  I[1, individual[i]] + I[2, opponent[i]];
   }
 }

model {
// Create vector of predicted values
  to_vector([B_0]) ~ normal(0, 1);

 // Random effects distribution
    to_vector(zI) ~ normal(0, 1);
    to_vector(sigma_I) ~ normal(0, 1);
    L ~ lkj_corr_cholesky(2);

 // Likelihood function
    for (i in 1:n_obs)
     z[i]~normal(e_z[i], sigma_e);
 }

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_phi;
real<lower=0> Sigma2_e;


real cov_2; //   int-resimpact
real cov_int_phi;


matrix[2, 2]  Omega_I;

Sigma2_intercept=sigma_I[1]^2;
Sigma2_phi=sigma_I[2]^2;
Sigma2_e=sigma_e^2;

Omega_I = L * L';
cov_2 = Omega_I[1,2]*sqrt(Sigma2_phi*Sigma2_intercept);
cov_int_phi = cov_2;


 }", file="RR.stan")

## data 800_400_2 ----------------------------
stan_func<-function(dfl_800_400_2_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_400_2_1x,
             pars = c("B_0", "Sigma2_intercept", "Sigma2_phi", "cov_int_phi", "cov_2", "Sigma2_e"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_400_2_1x <- do.call(rbind,parSapply(cl, dfl_800_400_2_1x, stan_func))
res_800_400_2_1x <- as.data.frame(res_800_400_2_1x)
res_800_400_2_1x$Parameter <- rep(c("B_0", "Sigma2_intercept", "Sigma2_phi", "cov_int_phi", "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_400_2_1x, "res_800_400_2_1x_M2.csv")


## data 800_200_4 ----------------------------
stan_func<-function(dfl_800_200_4_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_200_4_1x,
             pars = c("B_0", "Sigma2_intercept", "Sigma2_phi", "cov_int_phi", "cov_2", "Sigma2_e"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_200_4_1x <- do.call(rbind,parSapply(cl, dfl_800_200_4_1x, stan_func))
res_800_200_4_1x <- as.data.frame(res_800_200_4_1x)
res_800_200_4_1x$Parameter <- rep(c("B_0", "Sigma2_intercept", "Sigma2_phi", "cov_int_phi", "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_200_4_1x, "res_800_200_4_1x_M2.csv")


## data 800_100_8 ----------------------------
stan_func<-function(dfl_800_100_8_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_100_8_1x,
             pars = c("B_0", "Sigma2_intercept", "Sigma2_phi", "cov_int_phi", "cov_2", "Sigma2_e"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_100_8_1x <- do.call(rbind,parSapply(cl, dfl_800_100_8_1x, stan_func))
res_800_100_8_1x <- as.data.frame(res_800_100_8_1x)
res_800_100_8_1x$Parameter <- rep(c("B_0", "Sigma2_intercept", "Sigma2_phi", "cov_int_phi", "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_100_8_1x, "res_800_100_8_1x_M2.csv")


## data 800_50_16 ----------------------------
stan_func<-function(dfl_800_50_16_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_50_16_1x,
             pars = c("B_0", "Sigma2_intercept", "Sigma2_phi", "cov_int_phi", "cov_2", "Sigma2_e"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_50_16_1x <- do.call(rbind,parSapply(cl, dfl_800_50_16_1x, stan_func))
res_800_50_16_1x <- as.data.frame(res_800_50_16_1x)
res_800_50_16_1x$Parameter <- rep(c("B_0", "Sigma2_intercept", "Sigma2_phi", "cov_int_phi", "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_50_16_1x, "res_800_50_16_1x_M2.csv")

# M3 ---------------------------
write("data {
//Phenotype dataset

  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or aantal rows
   int<lower=0> n_ind; //number of individuals

  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent[n_obs];  //  Individual ID opponent repeated obs


  // Continuous variables
    real xj[n_obs];  // observed covariate (opponent trait) measured with error
     real z[n_obs];  // phenotypic observations
 }

 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
     real B_0; //intercept
     real psi; //slope


   // Random effects
   matrix[2,n_ind]         zI; //(intercepts and slopes, res_impact for each individual)
   vector<lower=0>[2]      sigma_I; // sd  intercepts, slopes, res_impact
   cholesky_factor_corr[2] L;  // factor to estimate covariance int-slopes
   real<lower=0> sigma_e;
 }

 transformed parameters{
    matrix[2,n_ind] I; //  Unscaled blups intercept and slope and res_impact for each individual
    real e_z[n_obs]; // predicted values for phenotype
    I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value

   for (i in 1:n_obs) {
     e_z[i]  = B_0  +  I[1, individual[i]] + psi*xj[i] + I[2, opponent[i]];
   }
 }

model {
// Create vector of predicted values
  B_0 ~ normal(0, 1);
  psi ~ normal(0, 1);

 // Random effects distribution
    to_vector(zI) ~ normal(0, 1);
    to_vector(sigma_I) ~ normal(0, 1);
    L ~ lkj_corr_cholesky(2);

 // Likelihood function
    for (i in 1:n_obs)
     z[i]~normal(e_z[i], sigma_e);

 }

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_epsilon;
real<lower=0> Sigma2_phi;
real<lower=0> Sigma2_e;


real cov_2; //   int-resimpact



matrix[2, 2]  Omega_I;

Sigma2_intercept=sigma_I[1]^2;
Sigma2_epsilon=sigma_I[2]^2;
Sigma2_phi=psi^2*variance(xj) + sigma_I[2]^2;
Sigma2_e=sigma_e^2;

Omega_I = L * L';
cov_2 = Omega_I[1,2]*sqrt(Sigma2_epsilon*Sigma2_intercept);


 }", file="RR.stan")

## data 800_400_2 ----------------------------
stan_func<-function(dfl_800_400_2_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_400_2_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_800_400_2_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                              "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                              "cov_2", "Sigma2_e", "lp__")

  }

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_400_2_1x <- do.call(rbind,parSapply(cl, dfl_800_400_2_1x_e, stan_func))
res_800_400_2_1x <- as.data.frame(res_800_400_2_1x)
res_800_400_2_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                             "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                             "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_400_2_1x, "res_800_400_2_1x_M3.csv")

## data 800_200_4 ----------------------------
stan_func<-function(dfl_800_200_4_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_200_4_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_800_200_4_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")
  
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_200_4_1x <- do.call(rbind,parSapply(cl, dfl_800_200_4_1x_e, stan_func))
res_800_200_4_1x <- as.data.frame(res_800_200_4_1x)
res_800_200_4_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_200_4_1x, "res_800_200_4_1x_M3.csv")




## data 800_100_8 ----------------------------
stan_func<-function(dfl_800_100_8_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_100_8_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_800_100_8_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")
  
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_100_8_1x <- do.call(rbind,parSapply(cl, dfl_800_100_8_1x_e, stan_func))
res_800_100_8_1x <- as.data.frame(res_800_100_8_1x)
res_800_100_8_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_100_8_1x, "res_800_100_8_1x_M3.csv")

## data 800_50_16 ----------------------------
stan_func<-function(dfl_800_50_16_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_50_16_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_800_50_16_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")
  
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_50_16_1x <- do.call(rbind,parSapply(cl, dfl_800_50_16_1x_e, stan_func))
res_800_50_16_1x <- as.data.frame(res_800_50_16_1x)
res_800_50_16_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_50_16_1x, "res_800_50_16_1x_M3.csv")



## M4 ----------------------------------
write("data {
//Phenotype dataset

  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or aantal rows
   int<lower=0> n_ind; //number of individuals

  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent[n_obs];  //  Individual ID opponent repeated obs


  // Continuous variables
    real xj[n_obs];  // observed covariate (opponent trait) measured with error
     real z[n_obs];  // phenotypic observations
 }

 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
     real B_x; //intercept covariate model
     real B_0; //intercept
     real psi; //slope


   // Random effects
   matrix[3,n_ind]         zI; //(intercepts and slopes, res_impact for each individual)
   vector<lower=0>[3]      sigma_I; // sd  intercepts, slopes, res_impact
   cholesky_factor_corr[3] L;  // factor to estimate covariance int-slopes
   real<lower=0> sigma_e;
   real<lower=0> sigma_ex;
 }

 transformed parameters{
    matrix[3,n_ind] I; //  Unscaled blups intercept and slope and res_impact for each individual
    real e_z[n_obs]; // predicted values for phenotype
    real e_x[n_obs]; // predicted values for covariate
    I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value

     for (i in 1:n_obs) {
     e_x[i] = B_x +  I[3, opponent[i]];
     }

   for (i in 1:n_obs) {
     e_z[i]  = B_0  +  I[1, individual[i]] + psi*I[3, opponent[i]] + I[2, opponent[i]];
   }
 }

model {
// Create vector of predicted values
  to_vector([B_0]) ~ normal(0, 1);
  to_vector([B_x]) ~ normal(0, 1);
  to_vector([psi]) ~ normal(0, 1);

 // Random effects distribution
    to_vector(zI) ~ normal(0, 1);
    to_vector(sigma_I) ~ normal(0, 1);
    L ~ lkj_corr_cholesky(3);

 // Likelihood function
    for (i in 1:n_obs)
     z[i]~normal(e_z[i], sigma_e);

     for (i in 1:n_obs)
     xj[i]~normal(e_x[i], sigma_ex);
 }

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_epsilon;
real<lower=0> Sigma2_x;
real<lower=0> Sigma2_phi;
real<lower=0> Sigma2_e;
real<lower=0> Sigma2_ex;

real cov_2; //   int-resimpact
real cov_4; //   int-x
real cov_6; //   x-resimpact
real cov_int_phi;


matrix[3, 3]  Omega_I;

Sigma2_intercept=sigma_I[1]^2;
Sigma2_epsilon=sigma_I[2]^2;
Sigma2_x=sigma_I[3]^2; 
Sigma2_e=sigma_e^2;
Sigma2_ex=sigma_ex^2;

Omega_I = L * L';
cov_2 = Omega_I[1,2]*sqrt(Sigma2_epsilon*Sigma2_intercept);
cov_4 = Omega_I[1,3]*sqrt(Sigma2_x*Sigma2_intercept);
cov_6 = Omega_I[2,3]*sqrt(Sigma2_x*Sigma2_epsilon);
cov_int_phi = cov_2 + psi*cov_4;
Sigma2_phi=sigma_I[2]^2 + psi^2*sigma_I[3]^2; // variance total impact, also depends on covariance between x and epsilon
 }", file="RR.stan")

## data 800_400_2 ----------------------------
stan_func<-function(dfl_800_400_2_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_400_2_1x,
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                       "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                       "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_400_2_1x <- do.call(rbind,parSapply(cl, dfl_800_400_2_1x, stan_func))
res_800_400_2_1x <- as.data.frame(res_800_400_2_1x)
res_800_400_2_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_800_400_2_1x, "res_800_400_2_1x_M4.csv")


## data 800_200_4 ----------------------------
stan_func<-function(dfl_800_200_4_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_200_4_1x,
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                       "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                       "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_200_4_1x <- do.call(rbind,parSapply(cl, dfl_800_200_4_1x, stan_func))
res_800_200_4_1x <- as.data.frame(res_800_200_4_1x)
res_800_200_4_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_800_200_4_1x, "res_800_200_4_1x_M4.csv")


## data 800_100_8 ----------------------------
stan_func<-function(dfl_800_100_8_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_100_8_1x,
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                       "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                       "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_100_8_1x <- do.call(rbind,parSapply(cl, dfl_800_100_8_1x, stan_func))
res_800_100_8_1x <- as.data.frame(res_800_100_8_1x)
res_800_100_8_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_800_100_8_1x, "res_800_100_8_1x_M4.csv")


## data 800_50_16 ----------------------------
stan_func<-function(dfl_800_50_16_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_50_16_1x,
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                       "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                       "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_50_16_1x <- do.call(rbind,parSapply(cl, dfl_800_50_16_1x, stan_func))
res_800_50_16_1x <- as.data.frame(res_800_50_16_1x)
res_800_50_16_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_800_50_16_1x, "res_800_50_16_1x_M4.csv")

## data 300 ----------------------------
stan_func<-function(dfl_300_75_4_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_300_75_4_1x,
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                      "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                      "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_300_75_4_1x <- do.call(rbind,parSapply(cl, dfl_300_75_4_1x, stan_func))
res_300_75_4_1x <- as.data.frame(res_300_75_4_1x)
res_300_75_4_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                         "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                         "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_300_75_4_1x, "res_300_75_4_1x_M4.csv")


## data 600 ----------------------------
stan_func<-function(dfl_600_150_4_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_600_150_4_1x,
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                      "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                      "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_600_150_4_1x <- do.call(rbind,parSapply(cl, dfl_600_150_4_1x, stan_func))
res_600_150_4_1x <- as.data.frame(res_600_150_4_1x)
res_600_150_4_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                             "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                             "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_600_150_4_1x, "res_600_150_4_1x_M4.csv")


## data 1200 ----------------------------
stan_func<-function(dfl_1200_300_4_1x){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_1200_300_4_1x,
             pars = c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                      "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                      "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex"),
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  dat
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_1200_300_4_1x <- do.call(rbind,parSapply(cl, dfl_1200_300_4_1x, stan_func))
res_1200_300_4_1x <- as.data.frame(res_1200_300_4_1x)
res_1200_300_4_1x$Parameter <- rep(c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                              "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                              "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex", "lp_"),iterations)
write.csv(res_1200_300_4_1x, "res_1200_300_4_1x_M4.csv")





## M5 --------------------------------
write("data {
//Phenotype dataset

  // Number of clusters (an integer)
   int<lower=0> n_obs; // number of observations for phenotypes or aantal rows
   int<lower=0> n_ind; //number of individuals

  // Clusters identifiers (an integer)
   int<lower=0> individual[n_obs];  //  Individual ID repeated obs
   int<lower=0> opponent[n_obs];  //  Individual ID opponent repeated obs

  // Continuous variables
    real xj[n_obs];  // observed covariate (opponent trait) measured with error
     real z[n_obs];  // phenotypic observations
 }

 parameters {
// Define parameters to estimate
  // Parameters for model of phenotype
     real B_0; //intercept
     real psi; //slope


   // Random effects
   matrix[3,n_ind]         zI; //(intercepts and slopes, res_impact for each individual)
   vector<lower=0>[3]      sigma_I; // sd  intercepts, slopes, res_impact
   cholesky_factor_corr[3] L;  // factor to estimate covariance int-slopes
   real<lower=0> sigma_e;
 }

 transformed parameters{
    matrix[3,n_ind] I; //  Unscaled blups intercept and slope and res_impact for each individual
    real e_z[n_obs]; // predicted values for phenotype
    I  = diag_pre_multiply(sigma_I, L) * zI; // get the unscaled value

   for (i in 1:n_obs) {
     e_z[i]  = B_0  +  I[1, individual[i]] + (psi+I[3, individual[i]])*xj[i] + I[2, opponent[i]];
   }
 }

model {
// Create vector of predicted values
  B_0 ~ normal(0, 1);
  psi ~ normal(0, 1);

 // Random effects distribution
    to_vector(zI) ~ normal(0, 1);
    to_vector(sigma_I) ~ normal(0, 1);
    L ~ lkj_corr_cholesky(3);

 // Likelihood function
    for (i in 1:n_obs)
     z[i]~normal(e_z[i], sigma_e);

 }

generated quantities{
real<lower=0> Sigma2_intercept;
real<lower=0> Sigma2_epsilon;
real<lower=0> Sigma2_psi;
real<lower=0> Sigma2_phi;
real<lower=0> Sigma2_e;


real cov_1; //   int-psi
real cov_2; //   int-resimpact
real cov_3; //   psi-resimpact


matrix[3,3]  Omega_I;


Sigma2_intercept=sigma_I[1]^2;
Sigma2_epsilon=sigma_I[2]^2;
Sigma2_psi=sigma_I[3]^2;
Sigma2_phi=psi^2*variance(xj) + sigma_I[2]^2;
Sigma2_e=sigma_e^2;


Omega_I = L * L';
cov_1 = Omega_I[1,3]*sqrt(Sigma2_psi*Sigma2_intercept);
cov_2 = Omega_I[1,2]*sqrt(Sigma2_epsilon*Sigma2_intercept);
cov_3 = Omega_I[2,3]*sqrt(Sigma2_psi*Sigma2_epsilon);


 }", file="RR.stan")

## data 800_400_2 ----------------------------
stan_func<-function(dfl_800_400_2_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_400_2_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_3000_1000_3_1x)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")
  
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_400_2_1x <- do.call(rbind,parSapply(cl, dfl_800_400_2_1x_e, stan_func))
res_800_400_2_1x <- as.data.frame(res_800_400_2_1x)
res_800_400_2_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_400_2_1x, "res_800_400_2_1x_M5.csv")

## data 800_200_4 ----------------------------
stan_func<-function(dfl_800_200_4_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_200_4_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_800_200_4_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")
  
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_200_4_1x <- do.call(rbind,parSapply(cl, dfl_800_200_4_1x_e, stan_func))
res_800_200_4_1x <- as.data.frame(res_800_200_4_1x)
res_800_200_4_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_200_4_1x, "res_800_200_4_1x_M5.csv")




## data 800_100_8 ----------------------------
stan_func<-function(dfl_800_100_8_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_100_8_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_800_100_8_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")
  
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_100_8_1x <- do.call(rbind,parSapply(cl, dfl_800_100_8_1x_e, stan_func))
res_800_100_8_1x <- as.data.frame(res_800_100_8_1x)
res_800_100_8_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_100_8_1x, "res_800_100_8_1x_M5.csv")

## data 800_50_16 ----------------------------
stan_func<-function(dfl_800_50_16_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_800_50_16_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_800_50_16_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")
  
}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_800_50_16_1x <- do.call(rbind,parSapply(cl, dfl_800_50_16_1x_e, stan_func))
res_800_50_16_1x <- as.data.frame(res_800_50_16_1x)
res_800_50_16_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                                    "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                                    "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_800_50_16_1x, "res_800_50_16_1x_M5.csv")




## data 300 ----------------------------
stan_func<-function(dfl_300_75_4_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_300_75_4_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_300_75_4_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi", "cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")

}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_300_75_4_1x <- do.call(rbind,parSapply(cl, dfl_300_75_4_1x_e, stan_func))
res_300_75_4_1x <- as.data.frame(res_300_75_4_1x)
res_300_75_4_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                             "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                             "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_300_75_4_1x, "res_300_75_4_1x_M5.csv")



## data 600 ----------------------------
stan_func<-function(dfl_600_150_4_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_600_150_4_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_600_150_4_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")

}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_600_150_4_1x <- do.call(rbind,parSapply(cl, dfl_600_150_4_1x_e, stan_func))
res_600_150_4_1x <- as.data.frame(res_600_150_4_1x)
res_600_150_4_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                             "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                             "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_600_150_4_1x, "res_600_150_4_1x_M5.csv")


## data 1200 ----------------------------
stan_func<-function(dfl_1200_300_4_1x_e){
  #Parameters monitored
  md <- stan("RR.stan", data = dfl_1200_300_4_1x_e,
             control=list(max_treedepth = 15), chains = 1, iter = 5000,
             warmup = 1000, thin = 1, cores = 1)
  dat <- as.data.frame(t(round(summary(md)$summary[,c(1,4,6,8,9,10)],3)))
  # calculate variance of opponent trait xj form simulated data
  df2 <- as.data.frame(dfl_1200_300_4_1x_e)
  df1 <- df2 %>% group_by(opponent) %>% slice_sample(n=1)
  var_xj <- var(df1$xj)
  # calculate covariance int-phi
  blups_int <- dat %>% select(starts_with("I[1,"))
  blups_int <- as.data.frame(t(blups_int))
  dat$cov_int_phi <- dat$cov_2 + dat$psi*(cor(blups_int$`50%`,df1$xj)*sqrt(dat$Sigma2_intercept*var_xj))
  dat <- select(dat, "B_0", "psi", "Sigma2_intercept",
                "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                "cov_2", "Sigma2_e", "lp__")

}

## Apply function to list of data sets
## Not that the mclapply function needs to be modified for windows.
res_1200_300_4_1x <- do.call(rbind,parSapply(cl, dfl_1200_300_4_1x_e, stan_func))
res_1200_300_4_1x <- as.data.frame(res_1200_300_4_1x)
res_1200_300_4_1x$Parameter <- rep(c("B_0", "psi", "Sigma2_intercept",
                             "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                             "cov_2", "Sigma2_e", "lp_"),iterations)
write.csv(res_1200_300_4_1x, "res_1200_300_4_1x_M5.csv")
