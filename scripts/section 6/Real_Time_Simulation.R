## Distributed as part of the supporting materials for the manuscript
## "Flexible marked spatio-temporal point processes with applications to event sequences from association football"
##
## Author: Santhosh Narayanan
## Date: 16 Oct 2022
## Licence: GPL 3
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! PROVIDED "AS IS"
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

# Section 6 Script 4
# Real time simulation

library(parallel)
library(coda)
library(data.table)
library(matrixStats)
library(parallel)
library(ggplot2)
library(pROC)

setwd('~/Documents/PhD/WorkingPaper/')
source('Section6/helpers.R')
load("matchData.RData")
source('Section6/Model_5_100_teams/data_file.R')

myfiles = lapply(paste0('Section6/Model_5_100_teams/sample_file_',c(1:3),'.csv'), read.csv, skip = 40)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) {
  df <- df[c(5:504), -c(1:7)]
  mcmc(df)}
))
colnames(samples) <-  dotfnames_to_sqrfnames(colnames(samples))

# Prepare params with thinning
M <- 30
G <- 40
params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 3),])
main.params = params

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

###### Time model samples
myfiles = lapply(paste0('Section6/Model_simple/sample_file_',c(1:3),'.csv'), read.csv, skip = 38)
time.samples <- do.call(rbind, lapply(myfiles,FUN = function(df) mcmc(df[c(1005:1504), 8:940])))
colnames(time.samples) <-  dotfnames_to_sqrfnames(colnames(time.samples))
time.samples = time.samples[seq(from = 1, to = nrow(time.samples), by = 3), 
                            (M*(M-1)+M+3):(M*(M-1)+M+32)] 
time.params = colMeans(time.samples)

###### Zone model samples
train_data[, p_mark := shift(mark, 1, type = 'lag'), by = 'id']
train_data[, p_zone := shift(zone_s, 1, type = 'lag'), by = 'id']
zone_data = na.omit(train_data[, .N, by = c('p_zone', 'p_mark', 'zone_s')])
zone_cast = dcast.data.table(zone_data, formula = p_zone + p_mark ~ zone_s, fill = 0, drop = F)
zone.params = as.data.frame(zone_cast)
# Add 1 to the counts for the non-informative dirichlet prior
zone.params = cbind(zone.params[,c(1,2)], (zone.params[,-c(1,2)]+1)/rowSums(zone.params[,-c(1,2)]+1))

rm(list = c('zone_cast', 'zone_data', 'myfiles'))

########### Shot event probs evolution
cl <- makeCluster(getOption("cl.cores", 15))
clusterExport(cl, list("game_simulator_teams", "name_check",'time.params', 'zone.params',
                       'B1' , 'indx_def', 'indx_mid', 'indx_atk', 'logSumExp'))

model.results.shots <- vector(mode = "list", length = 20)

for(i in 1:20){
  sim.id = seq(41, 80, 2)[i]
  sim.meta = unlist(meta[sim.id,-1])
  sim.data = data_season_1[id %in% c(sim.id, sim.id+1),]
  sim.data$time[sim.data$period == 2] <- sim.data$time[sim.data$period == 2] + max(sim.data$time[sim.data$period == 1]) + 10
  sim.data = sim.data[, .SD, .SDcols = c('time', 'zone_s' , 'mark', 'history')]
  sim.time = floor(max(sim.data$time)/30)
  
  set.seed(123)
  
  for(j in 16:sim.time){
    sim.history = data.frame(sim.data[sim.data$time < (j-1)*30,])
    if(nrow(sim.history) < 1){
      sim.history <- data.frame(sim.data[1,])
    }
    tab = table(sim.history$mark)
    attributes(tab)$class <- "vector" 
    history.counts <- name_check(tab, 30)
    
    clusterExport(cl, list("sim.history","j",'history.counts','sim.meta'))
    
    sim.results = parApply(cl, main.params, 1, FUN = function(sim.params){
      results = sapply(1:20, FUN = function(z){
        sims = game_simulator_teams(par.sample = sim.params, time.params = time.params, 
                                    zone.params =  zone.params, meta.g = sim.meta,
                                    history = sim.history, Tend = j*30, M = 30)
        tab = table(sims$mark)
        attributes(tab)$class <- "vector" 
        return(name_check(tab, 30) - history.counts)
      })
    })
    
    reshape.results <- data.frame(NA)
    for(k in 1:20){
      reshape.results <- cbind(reshape.results, sim.results[(30*(k-1) + 1):(30*k),])
    }
    reshape.results <- as.matrix(reshape.results[,-1])
    
    model.results.shots[[i]] <- append(model.results.shots[[i]], list(reshape.results) )
  }
}

stopCluster(cl)

event_id <- 5
df.plot <- data.frame()

for(i in 1:20){
  sim.id = seq(41, 80, 2)[i]
  sim.data = data_season_1[id %in% c(sim.id, sim.id+1),]
  sim.data$time[sim.data$period == 2] <- sim.data$time[sim.data$period == 2] + max(sim.data$time[sim.data$period == 1]) + 10
  sim.data = sim.data[, .SD, .SDcols = c('time', 'zone_s' , 'mark', 'history')]
  sim.time = floor(max(sim.data$time)/30)
  
  true_events <- rep(0, sim.time)
  true_times <- unique(ceiling(sim.data$time[sim.data$mark == event_id]/30))
  true_events[true_times] <- 1
  model_events_probs <- rep(0, sim.time - 15)
  for(j in 1:(sim.time - 15)){
    my.mat <- model.results.shots[[i]][[j]]
    my.mat[my.mat > 1] <- 1
    model_events_probs[j] <- rowMeans(my.mat)[event_id]
  }
  
  df.plot <- rbind(df.plot,
                   data.frame(x = 0:(sim.time - 16), 
                        Truth = true_events[16:(sim.time)], 
                        Model = model_events_probs,
                        MA_5 =  c(NA, frollmean(head(true_events, -1), 5, fill = NA))[16:(sim.time)],
                        MA_10 = c(NA, frollmean(head(true_events, -1), 10, fill = NA))[16:(sim.time)],
                        MA_15 = c(NA, frollmean(head(true_events, -1), 15, fill = NA))[16:(sim.time)])
  )
}  

ggrocs(rocs = list(Model = roc(df.plot$Truth, df.plot$Model, direction = '<'),
                   MA_5 = roc(df.plot$Truth, df.plot$MA_5, direction = '<'),
                   MA_10 = roc(df.plot$Truth, df.plot$MA_10, direction = '<'),
                   MA_15 = roc(df.plot$Truth, df.plot$MA_15, direction = '<')))

df.plot.melt <- melt.data.table(data.table(df.plot[1:102,-c(2, 4, 6)]), 1)
df.plot.melt$x  <- df.plot.melt$x + 15
df.plot.melt$xend <- df.plot.melt$x + 1
df.plot.melt$x <- df.plot.melt$x/2
df.plot.melt$xend <- df.plot.melt$xend/2
ggplot(data = df.plot.melt, aes(x = x, y = value, xend = xend, yend = value, colour = variable)) +
  geom_segment() +
  geom_vline(xintercept = (true_times[true_times > 15])/2 - 0.25 , linetype = 'dotted') +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = 'Match time in minutes split into intervals of 30 seconds', y = 'P(atleast 1 Home Shot)') + ylim(c(0, 0.5)) + xlim(c(7, 59)) +
  theme(legend.position = 'bottom', legend.title = element_blank())
