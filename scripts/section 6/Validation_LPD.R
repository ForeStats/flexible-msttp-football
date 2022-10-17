## Distributed as part of the supporting materials for the manuscript
## "Flexible marked spatio-temporal point processes with applications to event sequences from association football"
##
## Author: Santhosh Narayanan
## Date: 16 Oct 2022
## Licence: GPL 3
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! PROVIDED "AS IS"
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

# Section 6 Script 2
# Bayesian validation - log predictive density

library(data.table)
library(coda)
library(matrixStats)
library(MCMCpack)
library(parallel)

setwd('~/Documents/PhD/WorkingPaper/')
source('Section6/helpers.R')

# Extract training data
load("matchData.RData")
M <- 30
G <- 40
meta <- data.table(expand.grid(game = as.factor(1:380), period = as.factor(1:2)))
setkeyv(meta, c('game', 'period'))
meta$id <- 1:760
data_season_1 <- meta[data_season_1, on = c('game', 'period')]
train_data <- data_season_1[id <= G]

#
cl <- makeCluster(getOption("cl.cores", 15))
clusterEvalQ(cl, library("data.table"))
clusterEvalQ(cl, library("matrixStats"))

# Locations model
train_data[, p_mark := shift(mark, 1, type = 'lag'), by = 'id']
train_data[, p_zone := shift(zone_s, 1, type = 'lag'), by = 'id']
zone_data = na.omit(train_data[, .N, by = c('p_zone', 'p_mark', 'zone_s')])
transitions = as.data.frame(dcast.data.table(zone_data, formula = p_zone + p_mark ~ zone_s, fill = 0, drop = F))
zone.samples <- array(1/30, dim = c(100, 90, 3))

set.seed(123)
for(i in 1:90){
  zone.samples[ , i, ] = rdirichlet(100, unlist(transitions[i, -c(1,2)]) + 1)
}

val.res.zones.model = sapply(1:10, FUN = function(z) 
  model.likelihood.zones(params = zone.samples, period.data = data_season_1[id == 40+z], M = 30, cl = cl))

# Simple beta model
source('Section6/Model_simple/data_file.R')

myfiles = lapply(paste0('Section6/Model_simple/sample_file_',c(1:3),'.csv'), read.csv, skip = 38)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) mcmc(df[c(1005:1504), 8:940])))
colnames(samples) <-  dotfnames_to_sqrfnames(colnames(samples))

# Prepare params
set.seed(123)
post.params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 15),])

val.res.times.model = sapply(1:10, FUN = function(z) 
  model.likelihood.times(params = post.params, period.data = data_season_1[id == 40+z], M = 30, cl = cl))

val.res.single.beta = sapply(1:10, FUN = function(z) 
  model.likelihood.marks.simple(params = post.params, period.data = data_season_1[id == 40+z], M = 30, cl = cl))

## Merge times and zones
val.res.times.model <- mapply("+", val.res.zones.model, val.res.times.model, SIMPLIFY = FALSE)

# Total
sum(sapply(list(mapply("+", val.res.single.beta, val.res.times.model, SIMPLIFY = FALSE)), 
           FUN = function(y) do.call(rbind, lapply(y, FUN = function(z) sum(rowMeans(z))))))
# Marks only
sum(do.call(rbind, lapply(val.res.single.beta, FUN = function(z) sum(rowMeans(z)))))
# Times only
sum(do.call(rbind, lapply(val.res.times.model, FUN = function(z) sum(rowMeans(z)))))

# Vector beta model
source('Section6/Model_vector/data_file.R')

myfiles = lapply(paste0('Section6/Model_vector/sample_file_',c(1:3),'.csv'), read.csv, skip = 38)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) mcmc(df[c(1005:1504), 8:954])))
colnames(samples) <-  dotfnames_to_sqrfnames(colnames(samples))

# Prepare params
set.seed(123)
post.params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 15),])

val.res.vector.beta = sapply(1:10, FUN = function(z) 
  model.likelihood.marks.vector(params = post.params, period.data = data_season_1[id == 40+z], M = 30, cl = cl))

# Total
sum(sapply(list(mapply("+", val.res.vector.beta, val.res.times.model, SIMPLIFY = FALSE)), 
           FUN = function(y) do.call(rbind, lapply(y, FUN = function(z) sum(rowMeans(z))))))
# Marks only
sum(do.call(rbind, lapply(val.res.vector.beta, FUN = function(z) sum(rowMeans(z)))))

# Beta matrix models

# History 5 - Threshold 50
source('Section6/Model_5_50/data_file.R')

myfiles = lapply(paste0('Section6/Model_5_50/sample_file_',c(1:3),'.csv'), read.csv, skip = 38)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) mcmc(df[c(1005:1504), -c(1:7)])))
colnames(samples) <-  dotfnames_to_sqrfnames(colnames(samples))

# Prepare params
set.seed(123)
post.params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 15),])

indx_atk <- apply(indx_atk, 2, as.integer)
indx_mid <- apply(indx_mid, 2, as.integer)
indx_def <- apply(indx_def, 2, as.integer)

game.lik <- model.likelihood.marks(params = post.params, period.data = data_season_1[id == 41], M = 30, cl = cl)
plot(rowMeans(game.lik), pch = 19)

val.res.5.50 = sapply(1:10, FUN = function(z) 
  model.likelihood.marks(params = post.params, period.data = data_season_1[id == 40+z], M = 30, cl = cl))

# Total
sum(sapply(list(mapply("+", val.res.5.50, val.res.times.model, SIMPLIFY = FALSE)), 
           FUN = function(y) do.call(rbind, lapply(y, FUN = function(z) sum(rowMeans(z))))))
# Marks only
sum(do.call(rbind, lapply(val.res.5.50, FUN = function(z) sum(rowMeans(z)))))

# History 5 - Threshold 100
source('Section6/Model_5_100/data_file.R')

myfiles = lapply(paste0('Section6/Model_5_100/sample_file_',c(1:3),'.csv'), read.csv, skip = 38)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) mcmc(df[c(1005:1504), -c(1:7)])))
colnames(samples) <-  dotfnames_to_sqrfnames(colnames(samples))

# Prepare params
set.seed(123)
post.params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 15),])

indx_atk <- apply(indx_atk, 2, as.integer)
indx_mid <- apply(indx_mid, 2, as.integer)
indx_def <- apply(indx_def, 2, as.integer)

game.lik <- model.likelihood.marks(params = post.params, period.data = data_season_1[id == 41], M = 30, cl = cl)
plot(rowMeans(game.lik), pch = 19)

val.res.5.100 = sapply(1:10, FUN = function(z) 
  model.likelihood.marks(params = post.params, period.data = data_season_1[id == 40+z], M = 30, cl = cl))

# Total
sum(sapply(list(mapply("+", val.res.5.100, val.res.times.model, SIMPLIFY = FALSE)), 
           FUN = function(y) do.call(rbind, lapply(y, FUN = function(z) sum(rowMeans(z))))))
# Marks only
sum(do.call(rbind, lapply(val.res.5.100, FUN = function(z) sum(rowMeans(z)))))

# History 10 - Threshold 50

# Extract training data
load("matchData_10.RData")
M <- 30
G <- 40
meta <- data.table(expand.grid(game = as.factor(1:380), period = as.factor(1:2)))
setkeyv(meta, c('game', 'period'))
meta$id <- 1:760
data_season_1 <- meta[data_season_1, on = c('game', 'period')]
train_data <- data_season_1[id <= G]

source('Section6/Model_10_50/data_file.R')

myfiles = lapply(paste0('Section6/Model_10_50/sample_file_',c(1:3),'.csv'), read.csv, skip = 38)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) mcmc(df[c(1005:1504), -c(1:7)])))
colnames(samples) <-  dotfnames_to_sqrfnames(colnames(samples))

# Prepare params
set.seed(123)
post.params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 15),])

indx_atk <- apply(indx_atk, 2, as.integer)
indx_mid <- apply(indx_mid, 2, as.integer)
indx_def <- apply(indx_def, 2, as.integer)

game.lik <- model.likelihood.marks(params = post.params, period.data = data_season_1[id == 41], M = 30, cl = cl)
plot(rowMeans(game.lik), pch = 19)

val.res.10.50 = sapply(1:10, FUN = function(z) 
  model.likelihood.marks(params = post.params, period.data = data_season_1[id == 40+z], M = 30, cl = cl))

# Total
sum(sapply(list(mapply("+", val.res.10.50, val.res.times.model, SIMPLIFY = FALSE)), 
           FUN = function(y) do.call(rbind, lapply(y, FUN = function(z) sum(rowMeans(z))))))
# Marks only
sum(do.call(rbind, lapply(val.res.10.50, FUN = function(z) sum(rowMeans(z)))))

# History 10 - Threshold 100
source('Section6/Model_10_100/data_file.R')

myfiles = lapply(paste0('Section6/Model_10_100/sample_file_',c(1:3),'.csv'), read.csv, skip = 38)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) mcmc(df[c(1005:1504), -c(1:7)])))
colnames(samples) <-  dotfnames_to_sqrfnames(colnames(samples))

# Prepare params
set.seed(123)
post.params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 15),])

indx_atk <- apply(indx_atk, 2, as.integer)
indx_mid <- apply(indx_mid, 2, as.integer)
indx_def <- apply(indx_def, 2, as.integer)

game.lik <- model.likelihood.marks(params = post.params, period.data = data_season_1[id == 41], M = 30, cl = cl)
plot(rowMeans(game.lik), pch = 19)

val.res.10.100 = sapply(1:10, FUN = function(z)
  model.likelihood.marks(params = post.params, period.data = data_season_1[id == 40+z], M = 30, cl = cl))

# Total
sum(sapply(list(mapply("+", val.res.10.100, val.res.times.model, SIMPLIFY = FALSE)), 
           FUN = function(y) do.call(rbind, lapply(y, FUN = function(z) sum(rowMeans(z))))))
# Marks only
sum(do.call(rbind, lapply(val.res.10.100, FUN = function(z) sum(rowMeans(z)))))

## Teams models

# History 5 - Threshold 100 - Teams
source('Section6/Model_5_100_teams/data_file.R')

load("matchData.RData")
M <- 30
G <- 40
meta <- data.table(expand.grid(game = as.factor(1:380), period = as.factor(1:2)))
setkeyv(meta, c('game', 'period'))
meta$id <- 1:760
data_season_1 <- meta[data_season_1, on = c('game', 'period')]
train_data <- data_season_1[id <= G]

# team info
meta = meta[df_season_1, on = c('game')]
meta = meta[, .SD,.SDcols = c('id', 'period','h_id', 'a_id')]
setkey(meta, id)

myfiles = lapply(paste0('Section6/Model_5_100_teams/sample_file_',c(1:3),'.csv'), read.csv, skip = 40)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) mcmc(df[c(5:504), -c(1:7)])))
colnames(samples) <-  dotfnames_to_sqrfnames(colnames(samples))

# Prepare params
set.seed(123)
post.params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 15),])

indx_atk <- apply(indx_atk, 2, as.integer)
indx_mid <- apply(indx_mid, 2, as.integer)
indx_def <- apply(indx_def, 2, as.integer)

game.lik <- model.likelihood.marks.teams(params = post.params, period.data = data_season_1[id == 41], 
                                         meta.g = unlist(meta[41, 2:4]), M = 30, cl = cl)
plot(rowMeans(game.lik), pch = 19)

val.res.teams = sapply(1:10, FUN = function(z)
  model.likelihood.marks.teams(params = post.params, period.data = data_season_1[id == 40+z], 
                               meta.g = unlist(meta[40+z, 2:4]), M = 30, cl))

# Total
sum(sapply(list(mapply("+", val.res.teams, val.res.times.model, SIMPLIFY = FALSE)), 
           FUN = function(y) do.call(rbind, lapply(y, FUN = function(z) sum(rowMeans(z))))))
# Marks only
sum(do.call(rbind, lapply(val.res.teams, FUN = function(z) sum(rowMeans(z)))))

# History 5 - Threshold 100 - Teams - delta
source('Section6/Model_5_100_teams_delta/data_file.R')

load("matchData.RData")
M <- 30
G <- 40
meta <- data.table(expand.grid(game = as.factor(1:380), period = as.factor(1:2)))
setkeyv(meta, c('game', 'period'))
meta$id <- 1:760
data_season_1 <- meta[data_season_1, on = c('game', 'period')]
train_data <- data_season_1[id <= G]

# team info
meta = meta[df_season_1, on = c('game')]
meta = meta[, .SD,.SDcols = c('id', 'period','h_id', 'a_id')]
setkey(meta, id)

myfiles = lapply(paste0('Section6/Model_5_100_teams_delta/sample_file_',c(1:3),'.csv'), read.csv, skip = 40)
samples <- do.call(rbind, lapply(myfiles,FUN = function(df) mcmc(df[c(5:504), -c(1:7)])))
colnames(samples) <-  dotfnames_to_sqrfnames(colnames(samples))

# Prepare params
set.seed(123)
post.params = as.matrix(samples[seq(from = 1, to = nrow(samples), by = 15),])

indx_atk <- apply(indx_atk, 2, as.integer)
indx_mid <- apply(indx_mid, 2, as.integer)
indx_def <- apply(indx_def, 2, as.integer)

game.lik <- model.likelihood.marks.teams.delta(params = post.params, period.data = data_season_1[id == 41], 
                                         meta.g = unlist(meta[41, 2:4]), M = 30, cl = cl)
plot(rowMeans(game.lik), pch = 19)

val.res.teams.delta = sapply(1:10, FUN = function(z)
  model.likelihood.marks.teams.delta(params = post.params, period.data = data_season_1[id == 40+z], 
                               meta.g = unlist(meta[40+z, 2:4]), M = 30, cl = cl))

# Total
sum(sapply(list(mapply("+", val.res.teams.delta, val.res.times.model, SIMPLIFY = FALSE)), 
           FUN = function(y) do.call(rbind, lapply(y, FUN = function(z) sum(rowMeans(z))))))
# Marks only
sum(do.call(rbind, lapply(val.res.teams.delta, FUN = function(z) sum(rowMeans(z)))))

##  Baseline model ##

base.params <- matrix(0, 100, 90)
counts = as.numeric(table(train_data$mark, train_data$zone_s))
totaltime = sum(train_data[, max(time), by = id]$V1)

set.seed(123)
for(i in 1:90){
  base.params[, i] = rgamma(100, shape = counts[i] + 1, rate = totaltime + 1)
}

val.res.base = colMeans(sapply(1:10, FUN = function(z)
  base.likelihood(params = base.params, period.data = data_season_1[id == 40+z], M = 30)))

sum(val.res.base)

##  Smart Baseline model ##

base.params <- as.numeric(table(train_data$mark))
base.params <- base.params/sum(base.params)
train_data[, p_mark := shift(mark, 1, type = 'lag'), by = 'id']
transitions <- as.matrix.data.frame(table(train_data$p_mark, train_data$mark))
smartbase.params <- array(1/30, dim = c(100, 30, 30))

set.seed(123)
for(i in 1:30){
  smartbase.params[ , i, ] = rdirichlet(100, transitions[i, ] + 1)
}

val.res.smartbase = sapply(1:10, FUN = function(z)
  smartbase.likelihood.marks(params = smartbase.params, base.params = base.params, 
                             period.data = data_season_1[id == 40+z], M = 30))

# Total
sum(sapply(list(mapply("+", val.res.smartbase, val.res.times.model, SIMPLIFY = FALSE)), 
           FUN = function(y) do.call(rbind, lapply(y, FUN = function(z) sum(rowMeans(z))))))
# Marks only
sum(do.call(rbind, lapply(val.res.smartbase, FUN = function(z) sum(rowMeans(z)))))

## Compile results ##

all.res <- list(val.res.single.beta, val.res.vector.beta, val.res.5.50, 
                val.res.10.50, val.res.5.100, val.res.10.100, val.res.teams, val.res.teams.delta, val.res.smartbase)

df.res <- cbind(sapply(all.res, FUN = function(y) do.call(rbind, 
                                                          lapply( mapply("+", y, val.res.times.model, SIMPLIFY = FALSE), 
                                                                  FUN = function(z) mean(colSums(z))))), val.res.base)

colnames(df.res) <- c('Single_beta', 'Vector_beta', 'Hist_5_Thr_50', 'Hist_10_Thr_50', 'Hist_10_Thr_100', 
                      'Hist_5_Thr_100', 'Hist_5_Thr_100_teams', 'Hist_5_Thr_100_teams_constrained' , 'Smart_Baseline', 'Baseline')

rownames(df.res) <- paste0('game_', 1:10)
df.res <- rbind(df.res, colSums(df.res))
rownames(df.res)[11] <- 'Total'

library(kableExtra)
results.table <- data.frame(t(df.res))
results.table <- results.table[order(results.table$Total),]
results.table <- cbind(npar = c(90, 538, 538, 870, 902, 915, 1539, 1495, 988, 988), 
                       results.table)
kable(results.table[,c(1, 12)], booktabs = TRUE, caption = "My table", align = "c", 'latex', digits = 2)

kable(results.table, booktabs = TRUE, caption = "My table", align = "c", 'latex', digits = 1) %>%
  kableExtra::landscape()

stopCluster(cl)
