## Distributed as part of the supporting materials for the manuscript
## "Flexible marked spatio-temporal point processes with applications to event sequences from association football"
##
## Author: Santhosh Narayanan
## Date: 16 Oct 2022
## Licence: GPL 3
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! PROVIDED "AS IS"
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

# Section 6 Script 1
# MCMC Diagnostics and Parameter descriptions for the History = 5, N = 100, teams model

library(coda)
library(kableExtra)
library(corrplot)
library(bayesplot)
library(ggpubr)
library(gtools)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(viridis)
library(scales)

setwd('~/Documents/PhD/WorkingPaper/')
load("Section6/Model_5_100_teams/matchData.RData")
source('Section6/pairs_function.R')
source('Section6/Model_5_100_teams/data_file.R')
source('Section6/helpers.R')

myfiles <- lapply(paste0('Section6/Model_5_100_teams/sample_file_', c(1:3), '.csv'),
                 read.csv,
                 skip = 40)

samples <- do.call(rbind, lapply(
  myfiles,
  FUN = function(df) {
    df <- df[-c(1:4, 505:509),-c(1:7)]
    mcmc(df)
  }
))

fit <- mcmc.list(lapply(
  myfiles,
  FUN = function(df) {
    df <- df[-c(1:4, 505:509),-c(1:7)]
    mcmc(df)
  }
))

pars_oi <- colnames(samples)
pars_oi <- dotfnames_to_sqrfnames(pars_oi)
colnames(samples) <-  pars_oi

summary.fit <- data.frame(cbind(
  summary(fit)[[1]],
  gelman.diag(fit, multivariate = FALSE)[[1]],
  N_eff = effectiveSize(fit)
))

row.names(summary.fit) <- pars_oi

# Table 8
print.ind <- c(14, 65, 108, 141, 248, 275, 314, 365, 508, 541, 748, 775, 907, 908, 909, 991, 1038, 1039)
kable(
  round(summary.fit[print.ind, c(1, 2, 5, 7)], 2),
  booktabs = TRUE,
  caption = "My table",
  align = "c",
  'latex'
)

par.means <- colMeans(samples)
par.sd    <- apply(samples, 2, sd)

# Table 12
delta.mean <- t(matrix(par.means[(9*B1+1):(9*B1+(3*M))], nrow = 3))
delta.mean[delta.mean < 0.01] <- 0
delta.mean <- round(delta.mean, 2)
delta.sd <- t(matrix(par.sd[(9*B1+1):(9*B1+(3*M))], nrow = 3))
delta.sd[delta.sd < 0.01] <- 0
delta.sd <- round(delta.sd, 2)
delta <- matrix(paste0(delta.mean, " (", delta.sd, ")"), nrow = 3, ncol = 30, byrow = T)
delta <- t(delta)
delta[delta.mean == 0 | delta.sd == 0] <- "."
row.names(delta) <- mark_levels
kable(delta[1:15,], booktabs = TRUE, caption = "My table", align = "c", 'latex')
kable(delta[16:30,], booktabs = TRUE, caption = "My table", align = "c", 'latex')

gamma_par <- apply(samples, 1, FUN = function(par.sample){
  
  meta.g <- c(1, 11, 11)
  
  theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
  
  theta_def_raw  = par.sample[(3*B1+1):(5*B1)]
  theta_mid_raw  = par.sample[(5*B1+1):(7*B1)]
  theta_atk_raw  = par.sample[(7*B1+1):(9*B1)]
  
  for (i in 1:B1){
    theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i]
    if(indx_def[i, 1] < (M/2 + 1)){
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }else{
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }
    
    theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i]
    if(indx_mid[i, 1] < (M/2 + 1)){
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      } 
    }else{
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      }
    }
    
    theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i]
    if(indx_atk[i, 1] < (M/2 + 1)){
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }else{
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }
  }
  
  theta <- list(theta_zone_1, theta_zone_2, theta_zone_3)
  
  mu_atk  <- rbind(matrix(par.sample[(9*B1+3*M+2):(9*B1+3*M+1+19*(M-1))], nrow = 19), rep(0, M-1))
  
  mu.game <- c(mu_atk[meta.g[2], 1:(M/2)], mu_atk[meta.g[3], -c(1:(M/2))])
  
  gamma =  matrix(rep(mu.game, each = M), nrow = M, ncol = M-1)
  
  gamma_list <- rep(list(gamma),3)
  
  for(i in 1:3){
    
    gamma_iter <- theta[[i]] + gamma
    
    gamma_iter <- exp(gamma_iter)/(1+rowSums(exp(gamma_iter)))
    gamma_iter <- cbind(gamma_iter, 1 - rowSums(gamma_iter))
    
    gamma_list[[i]] <- gamma_iter
  }
  gamma_list[[2]]
})
gamma_post_means_zone_2 <- matrix(rowMeans(gamma_par), 30, 30)
gamma_post_means_zone_2[gamma_post_means_zone_2 < 0.01] <- 0
gamma_post_means_zone_2 <- round(gamma_post_means_zone_2, 2)
gamma_post_sd_zone_2 <- matrix(apply(gamma_par, 1, sd), 30, 30)
gamma_post_sd_zone_2[gamma_post_sd_zone_2 < 0.01] <- 0
gamma_post_sd_zone_2 <- round(gamma_post_sd_zone_2, 2)
gamma_par <- matrix(paste0(gamma_post_means_zone_2, " (", gamma_post_sd_zone_2, ")"), nrow = 30, ncol = 30)
gamma_par[gamma_post_means_zone_2 == 0 | gamma_post_sd_zone_2 == 0] <- "."
row.names(gamma_par) = colnames(gamma_par) = mark_levels
kbl(gamma_par, booktabs = TRUE, caption = "My table", align = "c", 'latex') %>% 
  kable_styling(latex_options = c("scale_down")) %>% row_spec(0, angle = 90)

# Figure 13
my.mat <- gamma_post_means_zone_2
colnames(my.mat) <- row.names(my.mat) <- mark_levels
plot.ind <- c(1,3,9)
plot.ind <- c(plot.ind, plot.ind + 15)
corrplot(my.mat[plot.ind, plot.ind],
         method="circle",
         is.corr=FALSE,
         type="full",
         tl.srt = 45,
         number.cex = 1,
         cl.lim=c(0,1), 
         addCoef.col = rgb(0,0,0, alpha = 0.6)
)

# Pass - Pass ability
gamma_par_team <- apply(samples, 1, FUN = function(par.sample){
  
  theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
  
  theta_def_raw  = par.sample[(3*B1+1):(5*B1)]
  theta_mid_raw  = par.sample[(5*B1+1):(7*B1)]
  theta_atk_raw  = par.sample[(7*B1+1):(9*B1)]
  
  for (i in 1:B1){
    theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i]
    if(indx_def[i, 1] < (M/2 + 1)){
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }else{
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }
    
    theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i]
    if(indx_mid[i, 1] < (M/2 + 1)){
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      } 
    }else{
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      }
    }
    
    theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i]
    if(indx_atk[i, 1] < (M/2 + 1)){
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }else{
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }
  }
  
  theta <- list(theta_zone_1, theta_zone_2, theta_zone_3)
  
  mu_atk  <- rbind(matrix(par.sample[(9*B1+3*M+2):(9*B1+3*M+1+19*(M-1))], nrow = 19), rep(0, M-1))
  
  difference_vec <- rep(0, 20)
  
  for(tm in 1:20){
    meta.g <- c(1, tm, tm)
    
    mu.game <- c(mu_atk[meta.g[2], 1:(M/2)], mu_atk[meta.g[3], -c(1:(M/2))])
    
    gamma =  matrix(rep(mu.game, each = M), nrow = M, ncol = M-1)
    
    gamma_iter <- theta[[2]] + gamma
    
    difference_vec[tm] <- gamma_iter[3, 3] - gamma_iter[18, 18]
  }
  difference_vec
})

team.names = sort(unique(df_season_1$home))

gamma_par_team <- t(gamma_par_team)
colnames(gamma_par_team) <- team.names
row.names(gamma_par_team) <- 1:nrow(gamma_par_team)
mu_melt <- melt(data.frame(gamma_par_team), variable.name =  'team', value.name = 'ability')
mu_melt$team <- sapply(mu_melt$team, FUN = function(z) gsub("\\.", " ", as.character(z)))
mu_melt$team <- factor(as.character(mu_melt$team), names(sort(colMeans(gamma_par_team))))

percents <- data.frame(team = team.names, pct = paste0(apply(gamma_par_team, 2, FUN = function(z) round(sum(z>0)/20, 0)), "%"))

ggplot(mu_melt,aes(y=team,x=ability,fill = ..x..))+
  geom_vline(xintercept = 0) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2) +
  scale_fill_viridis(name = "ability", option = "C", begin = 0.2) +
  xlim(c(-2.5,3.1)) +     # Use these limits for the Shot ability plot
  ylab('') + xlab('') +
  theme_bw(base_size=16) +
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.margin = margin(0.1, 0, -0.7, 0, "cm")) +
  theme(axis.ticks=element_blank()) + 
  annotate("text", x = 3, y = percents$team, label = percents$pct)

# Win - Pass ability
gamma_par_team <- apply(samples, 1, FUN = function(par.sample){
  
  theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
  
  theta_def_raw  = par.sample[(3*B1+1):(5*B1)]
  theta_mid_raw  = par.sample[(5*B1+1):(7*B1)]
  theta_atk_raw  = par.sample[(7*B1+1):(9*B1)]
  
  for (i in 1:B1){
    theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i]
    if(indx_def[i, 1] < (M/2 + 1)){
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }else{
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }
    
    theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i]
    if(indx_mid[i, 1] < (M/2 + 1)){
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      } 
    }else{
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      }
    }
    
    theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i]
    if(indx_atk[i, 1] < (M/2 + 1)){
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }else{
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }
  }
  
  theta <- list(theta_zone_1, theta_zone_2, theta_zone_3)
  
  mu_atk  <- rbind(matrix(par.sample[(9*B1+3*M+2):(9*B1+3*M+1+19*(M-1))], nrow = 19), rep(0, M-1))
  
  difference_vec <- rep(0, 20)
  
  for(tm in 1:20){
    meta.g <- c(1, tm, tm)
    
    mu.game <- c(mu_atk[meta.g[2], 1:(M/2)], mu_atk[meta.g[3], -c(1:(M/2))])
    
    gamma =  matrix(rep(mu.game, each = M), nrow = M, ncol = M-1)
    
    gamma_iter <- theta[[2]] + gamma
    
    difference_vec[tm] <- gamma_iter[1, 3] - gamma_iter[16, 18]
  }
  difference_vec
})

team.names = sort(unique(df_season_1$home))

gamma_par_team <- t(gamma_par_team)
colnames(gamma_par_team) <- team.names
row.names(gamma_par_team) <- 1:nrow(gamma_par_team)
mu_melt <- melt(data.frame(gamma_par_team), variable.name =  'team', value.name = 'ability')
mu_melt$team <- sapply(mu_melt$team, FUN = function(z) gsub("\\.", " ", as.character(z)))
mu_melt$team <- factor(as.character(mu_melt$team), names(sort(colMeans(gamma_par_team))))

percents <- data.frame(team = team.names, pct = paste0(apply(gamma_par_team, 2, FUN = function(z) round(sum(z>0)/20, 0)), "%"))

ggplot(mu_melt,aes(y=team,x=ability,fill = ..x..))+
  geom_vline(xintercept = 0) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2) +
  scale_fill_viridis(name = "ability", option = "C", begin = 0.2) +
  xlim(c(-2.5,3.6)) +     # Use these limits for the Shot ability plot
  ylab('') + xlab('') +
  theme_bw(base_size=16) +
  theme(panel.grid.minor=element_blank(), 
        panel.grid.major.x = element_blank(),
        plot.margin = margin(0.1, 0, -0.7, 0, "cm")) +
  theme(axis.ticks=element_blank()) + 
  annotate("text", x = 3.5, y = percents$team, label = percents$pct)

###
scaleFUN <- function(x) sprintf("%.1f", x)
print.ind <- c(14, 65, 108, 141, 248, 275, 314, 365, 508, 541, 748, 775, 907, 908, 909, 991, 1038, 1039)
mcmc_list_data <-  lapply(fit, FUN = function(df) {
  df <- df[, print.ind]
  colnames(df) <- c('beta[3%->%3*" | "*1]', 'beta[27%->%8*" | "*1]', 'beta[24%->%1*" | "*2]', 'beta[3%->%4*" | "*2]', 
                    'beta[3%->%5*" | "*3]', 'beta[3%->%10*" | "*3]', 'phi[3%->%3*" | "*1]', 'phi[27%->%8*" | "*1]', 
                    'phi[24%->%1*" | "*2]', 'phi[3%->%4*" | "*2]', 'phi[3%->%5*" | "*3]', 'phi[3%->%10*" | "*3]', 
                    'delta[3*" | "*1]', 'delta[3*" | "*2]', 'delta[3*" | "*3]', 'alpha', 'omega[9*", "*3]', 'omega[10*", "*3]')
  as.matrix(df)}
)

mcmc_trace(mcmc_list_data, pars = vars(1:18), facet_args = list(labeller=label_parsed, ncol = 3)) +
  xlab("Post-warmup iteration") + ylab("Parameter value") +
  scale_y_continuous(labels=scaleFUN)

print.ind <- c(14, 65, 108, 141, 248, 275, 314, 365, 508, 541, 748, 775, 907, 908, 909, 991, 1038, 1039)

#Posterior summary plot - SKIP
plot_data <- samples[, print.ind]
colnames(plot_data) <- c('beta[3%->%3*" | "*1]', 'beta[27%->%8*" | "*1]', 'beta[24%->%1*" | "*2]', 'beta[3%->%4*" | "*2]', 
                  'beta[3%->%5*" | "*3]', 'beta[3%->%10*" | "*3]', 'phi[3%->%3*" | "*1]', 'phi[27%->%8*" | "*1]', 
                  'phi[24%->%1*" | "*2]', 'phi[3%->%4*" | "*2]', 'phi[3%->%5*" | "*3]', 'phi[3%->%10*" | "*3]', 
                  'delta[3*" | "*1]', 'delta[3*" | "*2]', 'delta[3*" | "*3]', 'alpha', 'omega[9*", "*3]', 'omega[10*", "*3]')

plot.ind <- c(1:2, 7:8, 13:18)
pairsSN(plot_data[, plot.ind],
        panel =  function(...) points(..., pch = 21, col = '#011f4b', 
                                      add = T, bg = alpha("#011f4b", 0.75), cex = 0.75),
        diag.panel = panel.hist, cex.labels = 1.2, font.labels = 1)

######

prior.post.plot <- function(post_plot_data, stat_fun, stat_fun_args, ltitle, xlims, do = NULL){
  
  ggplot(post_plot_data) + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    stat_density(geom = 'line', position = 'identity', aes(x = value, colour = 'post', linetype = 'post'), n = 512) +
    #stat_function(fun = stat_fun, args = stat_fun_args, aes(colour = 'prior', linetype = 'prior')) +
    geom_line(data = data.frame(x = do$x, y = do$y), aes(x = x, y = y , colour = 'prior', linetype = 'prior')) +
    labs(title = ltitle, x = '', y ='') +
    scale_linetype_manual(values = c("post" = "solid",  "prior" = "dashed"),
                          name="distribution",
                          breaks=c("post", "prior"),
                          labels=c("posterior", "prior")) +
    scale_colour_manual(values=c("post" = "black", "prior" = "black"), 
                        name="distribution",
                        breaks=c("post", "prior"),
                        labels=c("posterior", "prior")) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_x_continuous(limits = xlims, breaks = pretty(xlims, 3), labels = pretty(xlims, 3)) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          plot.margin = margin(0.01, 0.01, 0.01, 0.01, "cm"))
}

# Prior overlay plots
print.ind <- c(14, 65, 108, 141, 248, 275, 314, 365, 508, 541, 748, 775, 907, 908, 909, 952, 953, 954, 1037, 1038, 1039, 1076, 1077, 1078)
# Beta

p1 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[1]]), 
                     stat_fun = dexp,
                     stat_fun_args = list(rate = 0.1),
                     ltitle = expression(beta[3%->%3*" | "*1]),
                     xlims = c(0, 5))

p2 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[2]]), 
                     stat_fun = dexp,
                     stat_fun_args = list(rate = 0.1),
                     ltitle = expression(beta[27%->%8*" | "*1]),
                     xlims = c(0, 5))

p3 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[3]]), 
                     stat_fun = dexp,
                     stat_fun_args = list(rate = 0.1),
                     ltitle = expression(beta[24%->%1*" | "*2]),
                     xlims = c(0, 5))

p4 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[4]]), 
                     stat_fun = dexp,
                     stat_fun_args = list(rate = 0.1),
                     ltitle = expression(beta[3%->%4*" | "*2]),
                     xlims = c(0, 5))

p5 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[5]]), 
                     stat_fun = dexp,
                     stat_fun_args = list(rate = 0.1),
                     ltitle = expression(beta[3%->%5*" | "*3]),
                     xlims = c(0, 5))

p6 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[6]]), 
                     stat_fun = dexp,
                     stat_fun_args = list(rate = 0.1),
                     ltitle = expression(beta[3%->%10*" | "*3]),
                     xlims = c(0, 5))

p7 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[7]]), 
                     stat_fun = dnorm,
                     stat_fun_args = list(mean = 0, sd = 10),
                     ltitle = expression(phi[3%->%3*" | "*1]),
                     xlims = c(-5, 5))

p8 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[8]]), 
                     stat_fun = dnorm,
                     stat_fun_args = list(mean = 0, sd = 10),
                     ltitle = expression(phi[27%->%8*" | "*1]),
                     xlims = c(-5, 5))

p9 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[9]]), 
                     stat_fun = dnorm,
                     stat_fun_args = list(mean = 0, sd = 10),
                     ltitle = expression(phi[24%->%1*" | "*2]),
                     xlims = c(-5, 5))

p10 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[10]]), 
                      stat_fun = dnorm,
                      stat_fun_args = list(mean = 0, sd = 10),
                     ltitle = expression(phi[3%->%4*" | "*2]),
                     xlims = c(-5, 5))

p11 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[11]]), 
                      stat_fun = dnorm,
                      stat_fun_args = list(mean = 0, sd = 10),
                     ltitle = expression(phi[3%->%5*" | "*3]),
                     xlims = c(-5, 5))

p12 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[12]]), 
                      stat_fun = dnorm,
                      stat_fun_args = list(mean = 0, sd = 10),
                     ltitle = expression(phi[3%->%10*" | "*3]),
                     xlims = c(-5, 5))

co <- rdirichlet(1000000, alpha = rep(1, 30))
do <- density(co[,1])

p13 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[13]]), 
                     stat_fun = ddirichlet,
                     stat_fun_args = list(min = 0, max = 1),
                     ltitle = expression(delta[3*" | "*1]),
                     xlims = c(0, 1),
                     do = do)

p14 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[14]]), 
                     stat_fun = dunif,
                     stat_fun_args = list(min = 0, max = 1),
                     ltitle = expression(delta[3*" | "*2]),
                     xlims = c(0, 1),
                     do = do)

p15 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[15]]), 
                      stat_fun = dunif,
                      stat_fun_args = list(min = 0, max = 1),
                      ltitle = expression(delta[3*" | "*3]),
                      xlims = c(0, 1),
                      do = do)

p16 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[16]]), 
                      stat_fun = ddirichlet,
                      stat_fun_args = list(min = 0, max = 1),
                      ltitle = expression(delta[5*" | "*1]),
                      xlims = c(0, 1),
                      do = do)

p17 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[17]]), 
                      stat_fun = dunif,
                      stat_fun_args = list(min = 0, max = 1),
                      ltitle = expression(delta[5*" | "*2]),
                      xlims = c(0, 1),
                      do = do)

p18 = prior.post.plot(post_plot_data = data.frame(value = samples[, print.ind[18]]), 
                      stat_fun = dunif,
                      stat_fun_args = list(min = 0, max = 1),
                      ltitle = expression(delta[5*" | "*3]),
                      xlims = c(0, 1),
                      do = do)

# Figure 8
figure = ggarrange(p1, p2, p3, p4, p5, p6, p7, p8,
                   p9, p10, p11, p12, p13, p14, p15, p16, p17, p18,
                   ncol=6, nrow=3, common.legend = TRUE, legend="top")

pdf("figures/postplot.pdf", width=10, height=8)
annotate_figure(figure,
                left = text_grob(label = "density", rot = 90),
                bottom = text_grob(label = "parameter value")
)
dev.off()

## High posterior density interval
library(hdrcde)
hdr(samples[,991])$hdr
hdr(exp(samples[,991]))$hdr

hdr(samples[,227])$hdr
hdr(samples[,253])$hdr

hdr(gamma_par[63,])$hdr
mean(gamma_par[63,])
hdr(gamma_par[528,])$hdr
mean(gamma_par[528,])
hdr(gamma_par[63,] - gamma_par[528,], prob = seq(50, 100, 1))$hdr
mean(gamma_par[63,] - gamma_par[528,])
hist(gamma_par[63,] - gamma_par[528,], n = 50)

## Figure 17 - Branching structure plot

# Prepare params with thinning
M <- 30
G <- 40
params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 20),])
main.params = params

# Extract training data
meta <- data.table(expand.grid(game = as.factor(1:380), period = as.factor(1:2)))
setkeyv(meta, c('game', 'period'))
meta$id <- 1:760
data_season_1 <- meta[data_season_1, on = c('game', 'period')]
train_data <- data_season_1[id <= G]
meta = meta[df_season_1, on = c('game')]
meta = meta[, .SD,.SDcols = c('id', 'period','h_id', 'a_id')]
setkey(meta, id)

indx_atk <- apply(indx_atk, 2, as.integer)
indx_mid <- apply(indx_mid, 2, as.integer)
indx_def <- apply(indx_def, 2, as.integer)

prob_array <- branch_str(main.params, train_data, meta, gid = 3, M = 30)

gid = 3  
events = train_data[id == gid]
n <- nrow(events)
b_str <- matrix(rowMeans(prob_array),n,n+1)[2:n, 1:n]
melted_bstr <- melt(b_str)
colnames(melted_bstr) <- c('Event', 'History', 'Probability')
melted_bstr <- melted_bstr[order(melted_bstr$Event),]
melted_bstr$History <- melted_bstr$History - 1
melted_bstr$Event_tag <- mark_levels[events$mark[melted_bstr$Event + 1]]
melted_bstr$History_tag[melted_bstr$History==0] <- 'Immigrant'
melted_bstr$History_tag[melted_bstr$History>0] <- mark_levels[events$mark[melted_bstr$History]]
melted_bstr$History[melted_bstr$History==0] <- melted_bstr$Event[melted_bstr$History==0]+1
melted_bstr <- melted_bstr[melted_bstr$Probability > 1e-2,]
melted_bstr$Lag <- -c(melted_bstr$Event - melted_bstr$History + 1)

melted_bstr <- melted_bstr[melted_bstr$Event >= 743,]
melted_bstr <- rbind(melted_bstr, list(742,	743, 1,	'Home_Pass_S',	'Immigrant', 0))
melted_bstr <- melted_bstr[order(melted_bstr$Event),]
melted_bstr <- melted_bstr[melted_bstr$Lag > -5,]
event.id = 783
ggplot(melted_bstr, aes(Event, Lag, fill = Probability)) + 
  geom_tile(colour = 'black', width=0.9, height=0.9, alpha = 0.8) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = 'bottom') +
  scale_fill_viridis_c(alpha = 0.8, begin = 0, end = 1,direction = -1, option = "A") + 
  scale_x_continuous(breaks = seq(min(melted_bstr$Event), max(melted_bstr$Event), 10),
                     labels = c(1, seq(10, 40, 10)),
                     position = "top") +
  scale_y_continuous(breaks=c(0:(-6)),
                     labels=c("Immigrant",
                              paste0('Offspring of ',c("i-1","i-2", "i-3", 'i-4', 'i-5', 'i-6')))) + 
  geom_label_repel(data=subset(melted_bstr, Event == event.id), aes(label = History_tag),
                   fill = 'blue', color = 'white', size = 2, segment.colour = 'black',
                   min.segment.length = 0.01) + 
  geom_label_repel(data=subset(melted_bstr, Event == event.id & Lag == -1), 
                   aes(x= Event, y = Lag + 2, label = Event_tag),
                   fill = 'red', color = 'white', size = 2, segment.colour = 'black',
                   min.segment.length = 0.01)  +
  labs(x = 'Event index (i)', y = 'Branching structure')
