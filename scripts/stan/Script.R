pmset -g
sudo pmset -a disablesleep 1
sudo pmset -a disablesleep 0

library(cmdstanr)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")

setwd("~/Documents/PhD/Stan/teamsHistory")
model.stan <- cmdstan_model("model.stan", cpp_options = list(stan_threads = TRUE))

source("data_file.R")
my_dat <- list(G = G,
               L = L,
               M = M,
               T = 20,
               B1 = B1,
               indx_def = indx_def,
               indx_mid = indx_mid,
               indx_atk = indx_atk,
               meta = meta,
               times = times,
               marks = marks)

initf <- function(chain_id = 1) {
  list(beta_def_raw  = rep(1, B1),
       beta_mid_raw  = rep(1, B1),
       beta_atk_raw  = rep(1, B1),
       theta_def_raw  = rep(0, 2*B1),
       theta_mid_raw  = rep(0, 2*B1),
       theta_atk_raw  = rep(0, 2*B1),
       delta = array(1/M, dim = c(3, M)),
       alpha = 0,
       mu_atk = array(0, dim = c(19, M-1)))
}
init_dat <- initf()

my_dat$grainsize <- 1
fit <- model.stan$sample(my_dat,
                         seed = 123,
                         refresh = 10,
                         init = list(init_dat),
                         output_dir = '~/Documents/PhD/Stan/teamsHistory',
                         chains = 1,
                         parallel_chains = 1,
                         threads_per_chain = 15,
                         iter_warmup = 500,
                         iter_sampling = 500,
                         save_warmup = TRUE)

save.image("workspace_teams.RData")

fit$save_object(file = "fit.RDS")
#fit2 <- readRDS("fit.RDS")

summary.df <- fit$summary()

inv_metric <- fit$inv_metric(matrix = F)[[1]]
step_size = 0.0568251

draws_array <- fit$draws()
draws_df <- as_draws_df(draws_array)
draws_df <- as.matrix(draws_df[,2:1543])
delta_mat <- matrix(draws_df[1, 901:990], nrow = 3, ncol = 30)
delta_mat[,30] <- 1 - rowSums(delta_mat[,-30])
  
init_dat <-  list(beta_def_raw  = draws_df[1, 1:100],
                  beta_mid_raw  = draws_df[1, 101:200],
                  beta_atk_raw  = draws_df[1, 201:300],
                  theta_def_raw  = draws_df[1, 301:500],
                  theta_mid_raw  = draws_df[1, 501:700],
                  theta_atk_raw  = draws_df[1, 701:900],
                  delta = delta_mat,
                  alpha = draws_df[1, 991],
                  mu_atk =  matrix(draws_df[1, 992:1542], nrow = 19, ncol = 29))

fit_resume <- model.stan$sample(my_dat,
                                seed = 1234,
                                refresh = 10,
                                init = rep(list(init_dat), 3),
                                output_dir = '~/Documents/PhD/Stan/teamsHistory',
                                chains = 3,
                                parallel_chains = 3,
                                threads_per_chain = 5,
                                iter_warmup = 0,
                                iter_sampling = 500,
                                adapt_engaged = FALSE,
                                step_size = step_size,
                                inv_metric = inv_metric)

fit_resume$save_object(file = "fit_resume.RDS")
summary.resume <- fit_resume$summary()
fit_resume$cmdstan_diagnose()
