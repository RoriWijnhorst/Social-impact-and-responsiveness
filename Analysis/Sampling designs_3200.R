library(rstan)
library(future)
library(future.apply)

# Compile model once
sm <- stan_model(file = "Models/I_R.stan")

# Parallel setup once
plan(multisession, workers = 50)

# Define stan_func once
stan_func <- function(data) {
  library(rstan)
  md <- sampling(sm, data = data,
                 pars = c("B_0", "B_x", "psi", "Sigma2_intercept", "Sigma2_psi", "Sigma2_x", "Sigma2_phi", 
                          "Sigma2_epsilon", "cov_int_psi", "cov_int_phi",
                          "cov_psi_phi","cor_1",  "cor_2", "cor_3", "cor_4", "cor_5",
                          "cor_6", "cov_1",  "cov_2", "cov_3", "cov_4", "cov_5",
                          "cov_6", "Sigma2_e", "Sigma2_ex"),
                 chains = 1, iter = 5000, warmup = 1000, thin = 1, cores = 1, save_warmup = FALSE)
  summ <- summary(md)$summary
  param_names <- rownames(summ)
  dat <- as.data.frame(round(summ[, c(1, 4, 6, 8, 9, 10)], 3))
  dat$Parameter <- param_names
  return(dat)
}

# --- 1. 
dfl <- readRDS("Data/dfl_3200_100_32_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_3200_100_32_1x_M1.csv", row.names = FALSE)

# --- 2. 
dfl <- readRDS("Data/dfl_3200_200_16_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_3200_200_16_1x_M1.csv", row.names = FALSE)

# --- 3. 
dfl <- readRDS("Data/dfl_3200_800_4_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_3200_800_4_1x_M1.csv", row.names = FALSE)

# --- 4. 
dfl <- readRDS("Data/dfl_3200_1600_2_1x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_3200_1600_2_1x_M1.csv", row.names = FALSE)

# --- 5. 
dfl <- readRDS("Data/dfl_3200_400_4_2x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_3200_400_4_2x_M1.csv", row.names = FALSE)

# --- 6. 
dfl <- readRDS("Data/dfl_3200_200_4_4x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_3200_200_4_4x_M1.csv", row.names = FALSE)

# --- 7. 
dfl <- readRDS("Data/dfl_3200_100_4_8x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_3200_100_4_8x_M1.csv", row.names = FALSE)

# --- 8.
dfl <- readRDS("Data/dfl_3200_200_8_2x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_3200_200_8_2x_M1.csv", row.names = FALSE)

# --- 9. 
dfl <- readRDS("Data/dfl_3200_200_2_8x.RDS")
output_list <- future_lapply(
  dfl, 
  stan_func, 
  future.globals = list(sm = sm),
  future.seed = TRUE
)
output <- do.call(rbind, output_list)
Sys.time()
write.csv(output, "Output/res_3200_200_2_8x_M1.csv", row.names = FALSE)