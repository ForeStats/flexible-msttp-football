## Distributed as part of the supporting materials for the manuscript
## "Flexible marked spatio-temporal point processes with applications to event sequences from association football"
##
## Author: Santhosh Narayanan
## Date: 16 Oct 2022
## Licence: GPL 3
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! PROVIDED "AS IS"
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

# Section 6 Script 3
# Team ability parameters

library(coda)
library(data.table)
library(matrixStats)
library(parallel)
library(ggplot2)
library(kableExtra)
library(tidyverse)
library(ggridges)
library(scales)
library(viridis)
library(reshape2)
library(xtable)
library(corrplot)

setwd('~/Documents/PhD/WorkingPaper/')
source('Section6/helpers.R')
load("Section6/Model_5_100_teams/matchData.RData")
source('Section6/Model_5_100_teams/data_file.R')

myfiles = lapply(paste0('Section6/Model_5_100_teams/sample_file_',c(1:3),'.csv'), read.csv, skip = 38)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) {
  df <- df[c(1005:1504), -c(1:7)]
  mcmc(df)}
))

pars_oi <- colnames(samples)
pars_oi <- dotfnames_to_sqrfnames(pars_oi)
colnames(samples) <-  pars_oi

# Conversion rates for Manchester City
gamma_par <- apply(samples, 1, FUN = function(par.sample){
  
  meta.g <- c(1, 11, 11)
    
  theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
  
  theta_def_raw  = par.sample[(3*B1+1+1):(3*B1+1+2*B1)]
  theta_mid_raw  = par.sample[(3*B1+1+2*B1+1):(3*B1+1+4*B1)]
  theta_atk_raw  = par.sample[(3*B1+1+4*B1+1):(3*B1+1+6*B1)]
  
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
  
  mu_atk  <- rbind(matrix(par.sample[(9*B1+1+3*(M-1)+1):(9*B1+1+22*(M-1))], nrow = 19), rep(0, M-1))
  omega   <- par.sample[(9*B1+1+22*(M-1)+1):(9*B1+1+23*(M-1))]
  
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

source("Section6/panel_function.R")

plot_data <- gamma_par
plot.ind <- c(1, 3, 9, 16, 18, 24)
plot.ind <- c(plot.ind, 
              plot.ind + 60,
              plot.ind + 240,
              plot.ind + 450,
              plot.ind + 510,
              plot.ind + 690)
panelSN(plot_data[, plot.ind],
        panel =  panel.hist,
        diag.panel = panel.hist)


gamma_post_means_zone_2 <- matrix(rowMeans(gamma_par), 30, 30)

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

# TBD Panel plot

######
# 
gamma_par_team <- apply(samples, 1, FUN = function(par.sample){
  
  theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
  
  theta_def_raw  = par.sample[(3*B1+1+1):(3*B1+1+2*B1)]
  theta_mid_raw  = par.sample[(3*B1+1+2*B1+1):(3*B1+1+4*B1)]
  theta_atk_raw  = par.sample[(3*B1+1+4*B1+1):(3*B1+1+6*B1)]
  
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
  
  mu_atk  <- rbind(matrix(par.sample[(9*B1+1+3*(M-1)+1):(9*B1+1+22*(M-1))], nrow = 19), rep(0, M-1))
  omega   <- par.sample[(9*B1+1+22*(M-1)+1):(9*B1+1+23*(M-1))]
  
  difference_vec <- rep(0, 20)
  
  for(tm in 1:20){
    meta.g <- c(1, tm, tm)
    
    mu.game <- c(mu_atk[meta.g[2], 1:(M/2)], mu_atk[meta.g[3], -c(1:(M/2))])
    
    gamma =  matrix(rep(mu.game, each = M), nrow = M, ncol = M-1)
    
    gamma_iter <- theta[[2]] + gamma
    
    gamma_iter <- exp(gamma_iter)/(1+rowSums(exp(gamma_iter)))
    gamma_iter <- cbind(gamma_iter, 1 - rowSums(gamma_iter))
    
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

ggplot(mu_melt,aes(y=team,x=ability,fill = ..x..))+
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2) +
  scale_fill_viridis(name = "ability", option = "C", begin = 0.2) +
  xlim(c(-0.3,0.3)) +     # Use these limits for the Shot ability plot
  theme_bw() + ylab('') + xlab('') +
  theme_bw(base_size=16) +
  theme(panel.grid.minor=element_blank(), plot.margin = margin(0.1, 0, -0.7, 0, "cm")) +
  theme(axis.ticks=element_blank())

gamma_par_team <- apply(samples, 1, FUN = function(par.sample){
  
  theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
  
  theta_def_raw  = par.sample[(3*B1+1+1):(3*B1+1+2*B1)]
  theta_mid_raw  = par.sample[(3*B1+1+2*B1+1):(3*B1+1+4*B1)]
  theta_atk_raw  = par.sample[(3*B1+1+4*B1+1):(3*B1+1+6*B1)]
  
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
  
  mu_atk  <- rbind(matrix(par.sample[(9*B1+1+3*(M-1)+1):(9*B1+1+22*(M-1))], nrow = 19), rep(0, M-1))
  omega   <- par.sample[(9*B1+1+22*(M-1)+1):(9*B1+1+23*(M-1))]
  
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

ggplot(mu_melt,aes(y=team,x=ability,fill = ..x..))+
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2) +
  scale_fill_viridis(name = "ability", option = "C", begin = 0.2) +
  xlim(c(-2.5,2.5)) +     # Use these limits for the Shot ability plot
  theme_bw() + ylab('') + xlab('') +
  theme_bw(base_size=16) +
  theme(panel.grid.minor=element_blank(), plot.margin = margin(0.1, 0, -0.7, 0, "cm")) +
  theme(axis.ticks=element_blank())

######

par.means = colMeans(samples)

mu_atk  <- rbind(matrix(par.means[(9*B1+1+3*(M-1)+1):(9*B1+1+22*(M-1))], nrow = 19), rep(0, M-1))

# Successful passing ability - Table 13
event_indx <- 3
team.names = sort(unique(df_season_1$home))
team.df = data.frame(team = team.names, Home_Pass_S = mu_atk[,event_indx], Away_Pass_S = mu_atk[,event_indx + 15])
team.df = team.df[order(team.df$Home_Pass_S + team.df$Away_Pass_S, decreasing = T),]
kable(team.df, booktabs = TRUE, caption = "My table", align = "c", 'latex', row.names = F, digits = 2)

# Team rankings figure 16
team.df =  data.frame(team = sort(unique(df_season_1$home)))
for(i in c( 3, 5, 10, 1, 7, 8, 11, 9)){
  event_indx <- i
  team.df = cbind(team.df, mu_atk[,event_indx] + mu_atk[,event_indx+15])
}
team.rankings <- data.frame(list(team = sort(unique(df_season_1$home)) , 
                                 21 - apply(team.df[,-1], 2, FUN = function(z) as.numeric(as.factor(z)))))
colnames(team.rankings) <- c('team',  'Pass', 'Shot', 'Goal', 'Win', 'Save', 'Clear', 'Foul', 'Lose')
team.rankings <- team.rankings[order(team.rankings$Pass),]

kable(team.rankings[,1:6], booktabs = TRUE, caption = "My table", align = "c", 'latex', row.names = F) %>% 
  add_header_above(c(" ", "Next zone"=6)) %>% 
  kable_styling(latex_options = "hold_position")

# Extract training data
G = 40
meta <- data.table(expand.grid(game = as.factor(1:380), period = as.factor(1:2)))
setkeyv(meta, c('game', 'period'))
meta$id <- 1:760
data_season_1 <- meta[data_season_1, on = c('game', 'period')]
train_data <- data_season_1[id <= G]
meta = meta[df_season_1, on = c('game')]
meta = meta[, .SD,.SDcols = c('id', 'period','h_id', 'a_id')]
setkey(meta, id)

#### shots in training per team
teamdata <- meta[train_data, on = c('id', 'period')]
homemarks = as.matrix.data.frame(table(teamdata$h_id, teamdata$mark))[,1:15]
awaymarks = as.matrix.data.frame(table(teamdata$a_id, teamdata$mark))[,-c(1:15)]
teammarks = homemarks + awaymarks
shotlist <- data.frame(team = team.names, shots = teammarks[,5])

homemarks = as.numeric(table(teamdata$h_id[teamdata$zone_s %in% c(3) & teamdata$mark %in% c(3, 4)]))
awaymarks = as.numeric(table(teamdata$a_id[teamdata$zone_s %in% c(1) & teamdata$mark %in% c(18, 19)]))
teamtouches = homemarks + awaymarks
shotlist$touches <- teamtouches
shotlist$shotpertouch <- shotlist$shots/shotlist$touches
shotlist <- shotlist[order(shotlist$shotpertouch, decreasing = T),]

kable(shotlist, booktabs = TRUE, caption = "My table", align = "c", 'latex', row.names = F, digits = 2)

# Figure 14
event_indx <- 3

mu_atk_hp  <- cbind(samples[,(9*B1+1+3*(M-1)+19*(event_indx-1)+1):(9*B1+1+3*(M-1)+19*(event_indx-1)+19)], 0)
mu_atk_ap  <- cbind(samples[,(9*B1+1+3*(M-1)+19*(event_indx+15-1)+1):(9*B1+1+3*(M-1)+19*(event_indx+15-1)+19)], 0)

# Use this for Figure 14 (a)
mu_atk <- mu_atk_hp

# Use this for Figure 14 (b)
mu_atk <- mu_atk_ap

# Figure 15
event_indx <- 5
mu_atk_hp  <- cbind(samples[,(9*B1+1+3*(M-1)+19*(event_indx-1)+1):(9*B1+1+3*(M-1)+19*(event_indx-1)+19)], 0)
mu_atk_ap  <- cbind(samples[,(9*B1+1+3*(M-1)+19*(event_indx+15-1)+1):(9*B1+1+3*(M-1)+19*(event_indx+15-1)+19)], 0)
mu_atk <- mu_atk_hp + mu_atk_ap

colnames(mu_atk) <- team.names
row.names(mu_atk) <- 1:nrow(mu_atk)
mu_melt <- melt(data.frame(mu_atk), variable.name =  'team', value.name = 'ability')
mu_melt$team <- sapply(mu_melt$team, FUN = function(z) gsub("\\.", " ", as.character(z)))
mu_melt$team <- factor(as.character(mu_melt$team), names(sort(colMeans(mu_atk))))

mu_melt$ability[mu_melt$team == 'West Ham United'] <- NA

ggplot(mu_melt,aes(y=team,x=ability,fill = ..x..))+
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01, quantile_lines = TRUE, quantiles = 2) +
  scale_fill_viridis(name = "ability", option = "C", begin = 0.2) +
  #xlim(c(-2.75,1.65)) +   # Use these limits for the Pass ability plot
  xlim(c(-0.8,1.6)) +     # Use these limits for the Shot ability plot
  theme_bw() + ylab('') + xlab('') +
  theme_bw(base_size=16) +
  theme(panel.grid.minor=element_blank(), plot.margin = margin(0.1, 0, -0.7, 0, "cm")) +
  theme(axis.ticks=element_blank()) +
  geom_point(aes(x = 0, y = 'West Ham United'))
