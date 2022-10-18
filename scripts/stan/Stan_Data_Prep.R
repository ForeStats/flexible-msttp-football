## Distributed as part of the supporting materials for the manuscript
## "Flexible marked spatio-temporal point processes with applications to event sequences from association football"
##
## Author: Santhosh Narayanan
## Date: 16 Oct 2022
## Licence: GPL 3
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! PROVIDED "AS IS"
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

## Script to prepare the input data and initialisations for HMC sampling via stan

library(rstan)
library(data.table)
rstan_options(auto_write = TRUE)

## Screening procedure based on Association Rules
load("Section6/Model_5_100_teams/matchData.RData")

G = 40
data_season_1 <- data_season_1[game %in% 1:G]
data_season_1 <- data.frame(data_season_1)

beta_matrix <- rep(list(matrix(0, 30, 30)),3)

for(i in 1:nrow(data_season_1)){
  if(data_season_1$history[i] > 0){
    
    zone_i <- data_season_1$zone_s[i]
    mark_i <- data_season_1$mark[i]
    mark_j <- unique(data_season_1$mark[(i - data_season_1$history[i]):(i - 1)])
    
    beta_matrix[[zone_i]][cbind(mark_j, mark_i)] <- 
      beta_matrix[[zone_i]][cbind(mark_j, mark_i)] + 1
  }
}

beta_compact <- rep(list(matrix(0, 15, 30)),3)

beta_compact[[1]] <-  beta_matrix[[1]][, 1:15] + rbind(beta_matrix[[3]][16:30, 16:30], beta_matrix[[3]][1:15, 16:30])
beta_compact[[2]] <-  beta_matrix[[2]][, 1:15] + rbind(beta_matrix[[2]][16:30, 16:30], beta_matrix[[2]][1:15, 16:30])
beta_compact[[3]] <-  beta_matrix[[3]][, 1:15] + rbind(beta_matrix[[1]][16:30, 16:30], beta_matrix[[1]][1:15, 16:30])

data_season_1 = data.table(data_season_1)

name.check <- function(vec,k){
  
  nam <- names(vec)
  add.k   <- which(!1:k %in% as.numeric(names(vec)))
  vec <- append(vec,rep(0,length(add.k)))
  names(vec) <- c(nam,add.k)
  vec <- vec[order(as.numeric(names(vec)))]
  names(vec) <- paste0('T',names(vec))
  
  return(vec)
}

# treatment zone 1

dummy = beta_compact[[1]][, -15]
zone_data = data_season_1[(mark %in% 1:14 & zone_s == 1) | (mark %in% 16:29 & zone_s == 3)]
total = nrow(zone_data)
xsupport = rowSums(dummy)/total
ysupport = name.check(table(zone_data$mark),29)
ysupport = c(ysupport[1:14] + ysupport[16:29])/total
xysupport = dummy/total

support = xysupport
sum(support > 0.0025, na.rm = T)
support[support < 0.0025] <- NA

lift = t(t(xysupport)/ysupport)/xsupport
sum(lift > 0.75, na.rm = T)
lift[lift < 0.75] <- NA

combined <- support*lift
sum(combined > 0, na.rm = T)
combined[combined < 0] <- NA
dummy[is.na(combined)] <- NA
sum(dummy > 0, na.rm = T)
min(dummy, na.rm = T)
beta_compact[[1]] <- dummy

# treatment zone 2

dummy = beta_compact[[2]][, -15]
zone_data = data_season_1[(mark %in% 1:14 & zone_s == 2) | (mark %in% 16:29 & zone_s == 2)]
total = nrow(zone_data)
xsupport = rowSums(dummy)/total
ysupport = name.check(table(zone_data$mark),29)
ysupport = c(ysupport[1:14] + ysupport[16:29])/total
xysupport = dummy/total

support = xysupport
sum(support > 0.0012, na.rm = T)
support[support < 0.0012] <- NA

lift = t(t(xysupport)/ysupport)/xsupport
sum(lift > 0.75, na.rm = T)
lift[lift < 0.75] <- NA

combined <- support*lift
sum(combined > 0, na.rm = T)
combined[combined < 0] <- NA
dummy[is.na(combined)] <- NA
sum(dummy > 0, na.rm = T)
min(dummy, na.rm = T)
beta_compact[[2]] <- dummy

# treatment zone 3

dummy = beta_compact[[3]][, -15]
zone_data = data_season_1[(mark %in% 1:14 & zone_s == 3) | (mark %in% 16:29 & zone_s == 1)]
total = nrow(zone_data)
xsupport = rowSums(dummy)/total
ysupport = name.check(table(zone_data$mark),29)
ysupport = c(ysupport[1:14] + ysupport[16:29])/total
xysupport = dummy/total

support = xysupport
sum(support > 0.00225, na.rm = T)
support[support < 0.00225] <- NA

lift = t(t(xysupport)/ysupport)/xsupport
sum(lift > 0.75, na.rm = T)
lift[lift < 0.75] <- NA

combined <- support*lift
sum(combined > 0, na.rm = T)
combined[combined < 0] <- NA
dummy[is.na(combined)] <- NA
sum(dummy > 0, na.rm = T)
min(dummy, na.rm = T)
beta_compact[[3]] <- dummy

rowSums(do.call(cbind, lapply(beta_compact, rowSums, na.rm = T)))

beta_melt <- lapply(beta_compact, FUN = melt, na.rm = T)
beta_melt <- lapply(beta_melt, FUN = function(z) z[z$value > 0,])
beta_melt <- lapply(beta_melt, FUN = function(z) z[order(z$Var2),-3])
beta_melt <- lapply(beta_melt, FUN = function(z) { 
  row.names(z) <- 1:nrow(z)
  return(z)
})

beta_def <- beta_melt[[1]]
beta_mid <- beta_melt[[2]]
beta_atk <- beta_melt[[3]]
B1 = nrow(beta_def)
B2 = nrow(beta_def[beta_def$Var2 < 15,])
B1 == B2
## END Screening

#Data set up
#Load data
load("Section6/Model_5_100_teams/matchData.RData")
G = 40

meta <- data.table(expand.grid(game = as.factor(1:380), period = as.factor(1:2)))
setkeyv(meta, c('game','period'))
meta$id <- 1:760
event_data <- meta[data_season_1, on = c('game', 'period')]

M = 30
g_len = array(0, dim = G)

for(i in 1:(G)){
  subset = event_data[id == i]
  subset[, length := seq(.N), by = id]
  g_len[i] = subset[, max(length)]
}

L = max(g_len)
times = array(999, dim=c(G, L))
marks = array(999, dim=c(G, 3*L))

for(i in 1:G){
  subset = event_data[id == i]
  times[i,1:nrow(subset)] = subset$time
  marks[i,1:nrow(subset)] = subset$history
  marks[i,(L+1):(L+nrow(subset))] = subset$zone_s
  marks[i,(2*L+1):(2*L+nrow(subset))] = subset$mark
}

metaData <- meta[df_season_1, on = 'game']
metaData$period <- as.numeric(metaData$period)
metaData$h_id <- as.numeric(metaData$h_id)
metaData$a_id <- as.numeric(metaData$a_id)
cols <- c("id", 'period','h_id',"a_id")
metaData <- metaData[,.SD,.SDcols = cols]
setkey(metaData, id)

my_dat <- list(G = G,
               L = L,
               M = M,
               T = 20,
               B1 = B1,
               indx_def = as.matrix(beta_def),
               indx_mid = as.matrix(beta_mid),
               indx_atk = as.matrix(beta_atk),
               meta = as.matrix(metaData[1:G, 2:4]),
               times = times,
               marks = marks)

stan_rdump(ls(my_dat), "data_file.R", envir = list2env(my_dat))

initf <- function(chain_id = 1) {
  list(beta_def_raw  = rep(1, B1),
       beta_mid_raw  = rep(1, B1),
       beta_atk_raw  = rep(1, B1),
       alpha = 0,
       theta_def_raw_sqd  = rep(0, 2*B1),
       theta_mid_raw_sqd  = rep(0, 2*B1),
       theta_atk_raw_sqd  = rep(0, 2*B1),
       delta_raw = array(0, dim = c(M-1, 3)),
       mu_atk_raw = array(0, dim = c(19, M-1)),
       omega_raw = rep(0, M-1),
       sigma_gamma = runif(1, 0, 5))
}
init <- initf()

stan_rdump(ls(init), "init_file.R", envir = list2env(init))
