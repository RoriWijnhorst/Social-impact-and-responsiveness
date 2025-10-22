library(rstan)
library(future)
library(future.apply)

# Compile model once
sm <- stan_model(file = "Models/VP.stan")

# Parallel setup once
plan(multisession, workers = 56)

# Define stan_func once
stan_func <- function(data) {
  library(rstan)
  md <- sampling(sm, data = data,
                 pars = c("B_0", "Sigma2_intercept", "Sigma2_phi", "cor_2", "cov_int_phi",
              "cov_2", "Sigma2_e"),
                 chains = 1, iter = 5000, warmup = 1000, thin = 1, cores = 1, save_warmup = FALSE)
  summ <- summary(md)$summary
  param_names <- rownames(summ)
  dat <- as.data.frame(round(summ[, c(1, 4, 6, 8, 9, 10)], 3))
  dat$Parameter <- param_names
  return(dat)
}


# --- 1. 
dfl <- readRDS("Data/dfl_800_50_16_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_50_16_1x_M2.csv", row.names = FALSE)

# --- 2. 
dfl <- readRDS("Data/dfl_800_100_8_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_100_8_1x_M2.csv", row.names = FALSE)

# --- 3. 
dfl <- readRDS("Data/dfl_800_200_4_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_200_4_1x_M2.csv", row.names = FALSE)

# --- 4. 
dfl <- readRDS("Data/dfl_800_400_2_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_400_2_1x_M2.csv", row.names = FALSE)


# Compile model once
sm <- stan_model(file = "Models/Trait.stan")

# Parallel setup once
plan(multisession, workers = 56)

# Define stan_func once
stan_func <- function(data) {
  library(rstan)
  md <- sampling(sm, data = data,
                 pars = c("B_0", "psi", "Sigma2_intercept",
                          "Sigma2_epsilon", "Sigma2_phi","Sigma2_x", "cov_int_phi",
                          "cov_2", "Sigma2_e"),
                 chains = 1, iter = 5000, warmup = 1000, thin = 1, cores = 1, save_warmup = FALSE)
  summ <- summary(md)$summary
  param_names <- rownames(summ)
  dat <- as.data.frame(round(summ[, c(1, 4, 6, 8, 9, 10)], 3))
  dat$Parameter <- param_names
  return(dat)
}

# --- 1. 
dfl <- readRDS("Data/dfl_800_50_16_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_50_16_1x_M3.csv", row.names = FALSE)

# --- 2. 
dfl <- readRDS("Data/dfl_800_100_8_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_100_8_1x_M3.csv", row.names = FALSE)

# --- 3. 
dfl <- readRDS("Data/dfl_800_200_4_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_200_4_1x_M3.csv", row.names = FALSE)

# --- 4. 
dfl <- readRDS("Data/dfl_800_400_2_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_400_2_1x_M3.csv", row.names = FALSE)



# Compile model once
sm <- stan_model(file = "Models/Trait_EIV.stan")

# Parallel setup once
plan(multisession, workers = 56)

# Define stan_func once
stan_func <- function(data) {
  library(rstan)
  md <- sampling(sm, data = data,
                 pars = c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_x",
                          "Sigma2_epsilon", "Sigma2_phi","cov_int_phi",
                          "cov_2", "cov_4", "cov_6","Sigma2_e", "Sigma2_ex"),
                 chains = 1, iter = 5000, warmup = 1000, thin = 1, cores = 1, save_warmup = FALSE)
  summ <- summary(md)$summary
  param_names <- rownames(summ)
  dat <- as.data.frame(round(summ[, c(1, 4, 6, 8, 9, 10)], 3))
  dat$Parameter <- param_names
  return(dat)
}

# --- 1. 
dfl <- readRDS("Data/dfl_800_50_16_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_50_16_1x_M4.csv", row.names = FALSE)

# --- 2. 
dfl <- readRDS("Data/dfl_800_100_8_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_100_8_1x_M4.csv", row.names = FALSE)

# --- 3. 
dfl <- readRDS("Data/dfl_800_200_4_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_200_4_1x_M4.csv", row.names = FALSE)

# --- 4. 
dfl <- readRDS("Data/dfl_800_400_2_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_400_2_1x_M4.csv", row.names = FALSE)


# Compile model once
sm <- stan_model(file = "Models/Trait_RS.stan")

# Parallel setup once
plan(multisession, workers = 56)

# Define stan_func once
stan_func <- function(data) {
  library(rstan)
  md <- sampling(sm, data = data,
                 pars = c("B_0", "psi", "Sigma2_intercept",
                          "Sigma2_epsilon", "Sigma2_phi", "Sigma2_x", "cov_int_phi",
                          "cov_2", "Sigma2_e"),
                 chains = 1, iter = 5000, warmup = 1000, thin = 1, cores = 1, save_warmup = FALSE)
  summ <- summary(md)$summary
  param_names <- rownames(summ)
  dat <- as.data.frame(round(summ[, c(1, 4, 6, 8, 9, 10)], 3))
  dat$Parameter <- param_names
  return(dat)
}


# --- 1. 
dfl <- readRDS("Data/dfl_800_50_16_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_50_16_1x_M5.csv", row.names = FALSE)

# --- 2. 
dfl <- readRDS("Data/dfl_800_100_8_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_100_8_1x_M5.csv", row.names = FALSE)

# --- 3. 
dfl <- readRDS("Data/dfl_800_200_4_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_200_4_1x_M5.csv", row.names = FALSE)

# --- 4. 
dfl <- readRDS("Data/dfl_800_400_2_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_800_400_2_1x_M5.csv", row.names = FALSE)