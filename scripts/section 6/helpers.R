## Distributed as part of the supporting materials for the manuscript
## "Flexible marked spatio-temporal point processes with applications to event sequences from association football"
##
## Author: Santhosh Narayanan
## Date: 16 Oct 2022
## Licence: GPL 3
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! PROVIDED "AS IS"
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

# Section 6 helpers

dotfnames_to_sqrfnames <- function(fnames) {
  fnames <- sapply(fnames,
                   function(i) {
                     if (!grepl("\\.", i)) return(i)
                     i <- sub("\\.", "[", i)
                     i <- sub("\\s*$", "]", i)
                     i }, USE.NAMES = FALSE)
  gsub("\\.\\s*", ",", fnames)
}

branch_str <- function(params, train_data, meta, gid, M){
  
  meta.g <- unlist(meta[gid, 2:4])
    
  prob_array <- apply(params, 1, FUN = function(par.sample){
    
    beta_zone_1 <-  beta_zone_2 <-  beta_zone_3 <- matrix(100, M, M)
    theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
    
    beta_def_raw  = par.sample[1:B1]
    beta_mid_raw  = par.sample[(B1+1):(2*B1)]
    beta_atk_raw  = par.sample[(2*B1+1):(3*B1)]
    theta_def_raw  = par.sample[(3*B1+1):(3*B1+2*B1)]
    theta_mid_raw  = par.sample[(3*B1+2*B1+1):(3*B1+4*B1)]
    theta_atk_raw  = par.sample[(3*B1+4*B1+1):(3*B1+6*B1)]
    
    for (i in 1:B1){
      beta_zone_1[indx_def[i, 1], indx_def[i, 2]] = beta_def_raw[i]
      theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i]
      if(indx_def[i, 1] < (M/2 + 1)){
        beta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
        if(indx_def[i, 2] < (M/2)){
          theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
        }
      }else{
        beta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
        if(indx_def[i, 2] < (M/2)){
          theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
        }
      }
      
      beta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = beta_mid_raw[i]
      theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i]
      if(indx_mid[i, 1] < (M/2 + 1)){
        beta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
        if(indx_mid[i, 2] < (M/2)){
          theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
        } 
      }else{
        beta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
        if(indx_mid[i, 2] < (M/2)){
          theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
        }
      }
      
      beta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = beta_atk_raw[i]
      theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i]
      if(indx_atk[i, 1] < (M/2 + 1)){
        beta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
        if(indx_atk[i, 2] < (M/2)){
          theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
        }
      }else{
        beta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
        if(indx_atk[i, 2] < (M/2)){
          theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
        }
      }
    }
    
    beta_list <- list(beta_zone_1, beta_zone_2, beta_zone_3)
    theta <- list(theta_zone_1, theta_zone_2, theta_zone_3)
    
    delta = matrix(par.sample[(9*B1+1):(9*B1+(3*M))], nrow = 3)
    alpha = par.sample[(9*B1+(3*M)+1)]
    
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
    
    events = train_data[id == gid]
    n <- nrow(events)
    
    prob_ij = cbind(1, matrix(0, nrow = n, ncol = n))
    
    for(i in 2:n){
      
      ti = events$time[i]
      mi = events$mark[i]
      zi = events$zone_s[i]
      hi = events$history[i]
      
      tj_vec = events$time[(i-hi):(i-1)]
      mj_vec = events$mark[(i-hi):(i-1)]
      
      gamma = gamma_list[[zi]]
      beta = beta_list[[zi]]
      
      g_ij_vec = gamma[mj_vec, mi]*exp(alpha-beta[mj_vec,mi]*(ti - tj_vec))
      mu_i = delta[zi,mi]
      g_ij_vec = c(mu_i, rep(0, i-hi-1), g_ij_vec, rep(0, n - i + 1))
      
      prob_ij[i,] = g_ij_vec/sum(g_ij_vec)
    }
    prob_ij
  })
}

name_check <- function(vec, k){
  
  nam <- names(vec)
  add.k   <- which(!1:k %in% as.numeric(nam))
  vec <- append(vec,rep(0,length(add.k)))
  names(vec) <- c(nam, add.k)
  vec <- vec[order(as.numeric(names(vec)))]
  names(vec) <- NULL
  
  return(vec)
}

game_simulator <- function(par.sample, time.params, zone.params, history, Tend, M){
  
  beta_zone_1 <-  beta_zone_2 <-  beta_zone_3 <- matrix(100, M, M)
  theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
  
  beta_def_raw  = par.sample[1:B1]
  beta_mid_raw  = par.sample[(B1+1):(2*B1)]
  beta_atk_raw  = par.sample[(2*B1+1):(3*B1)]
  theta_def_raw  = par.sample[(3*B1+1):(5*B1)]
  theta_mid_raw  = par.sample[(5*B1+1):(7*B1)]
  theta_atk_raw  = par.sample[(7*B1+1):(9*B1)]
  
  for (i in 1:B1){
    beta_zone_1[indx_def[i, 1], indx_def[i, 2]] = beta_def_raw[i]
    theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i]
    if(indx_def[i, 1] < (M/2 + 1)){
      beta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }else{
      beta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }
    
    beta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = beta_mid_raw[i]
    theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i]
    if(indx_mid[i, 1] < (M/2 + 1)){
      beta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      } 
    }else{
      beta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      }
    }
    
    beta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = beta_atk_raw[i]
    theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i]
    if(indx_atk[i, 1] < (M/2 + 1)){
      beta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }else{
      beta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }
  }
  
  beta_list <- list(beta_zone_1, beta_zone_2, beta_zone_3)
  theta <- list(theta_zone_1, theta_zone_2, theta_zone_3)
  
  delta = matrix(par.sample[(9*B1+1):(9*B1+(3*M))], nrow = 3)
  alpha = par.sample[(9*B1+(3*M)+1)]
  
  gamma_list <- rep(list(matrix(0, nrow = M, ncol = M-1)),3)
  
  for(i in 1:3){
    gamma_iter <- theta[[i]]
    
    gamma_iter <- exp(gamma_iter)/(1+rowSums(exp(gamma_iter)))
    gamma_iter <- cbind(gamma_iter, 1 - rowSums(gamma_iter))
    
    gamma_list[[i]] <- gamma_iter
  }
  
  shape = time.params[1:15]
  rate = time.params[16:30]
  
  # Initialization
  mark.history <- tail(history, 10)
  mark.n <- nrow(mark.history)
  
  n <- nrow(history)
  last_mark  <- tail(history$mark, 1)
  last_zone  <- tail(history$zone, 1)
  last_score  <- tail(history$score, 1)
  last_history <- tail(history$history, 1)
  
  next_zone <- sample(1:3, 1, prob = 
                        unlist(zone.params[zone.params$p_zone == last_zone & zone.params$p_mark == last_mark, 3:5]))
  
  next_time <- tail(history$time, 1) + rgamma(1, rep(shape, 2)[last_mark], rep(rate, 2)[last_mark])
  
  next_history = ifelse(last_mark %in% c(10:15, 25:30), 1, min(last_history+1, 5))
  
  gamma_mat = gamma_list[[next_zone]]
  beta_mat = beta_list[[next_zone]]
  
  loglik <- rep(0, M)
  
  for(k in 1:M){
    loglik[k] = logSumExp(c(log(delta[next_zone, k]),
                            alpha - beta_mat[tail(history$mark, next_history), k] * (next_time - tail(history$time, next_history) ) + 
                              log(gamma_mat[tail(history$mark, next_history), k]) ) )
  }
  mark_prob = exp(loglik - logSumExp(loglik))
  
  if(sum(mark_prob) == 0) mark_prob = gamma_mat[history$mark[n], ]
  
  next_mark = sample(1:M, 1, prob = mark_prob)
  
  history = rbind(history, list(next_time, next_zone, next_mark, next_history))
  n <- n + 1
  
  while(next_time < Tend){
    
    last_mark  <- tail(history$mark, 1)
    last_zone  <- tail(history$zone, 1)
    last_score  <- tail(history$score, 1)
    last_history <- tail(history$history, 1)
    
    next_zone <- sample(1:3, 1, prob = 
                          unlist(zone.params[zone.params$p_zone == last_zone & zone.params$p_mark == last_mark, 3:5]))
    
    next_time <- tail(history$time, 1) + rgamma(1, rep(shape, 2)[last_mark], rep(rate, 2)[last_mark])
    
    next_history = ifelse(last_mark %in% c(10:15, 25:30), 1, min(last_history+1, 5))
    
    gamma_mat = gamma_list[[next_zone]]
    beta_mat = beta_list[[next_zone]]
    
    loglik <- rep(0, M)
    
    for(k in 1:M){
      loglik[k] = logSumExp(c(log(delta[next_zone, k]),
                              alpha - beta_mat[tail(history$mark, next_history), k] * (next_time - tail(history$time, next_history) ) + 
                                log(gamma_mat[tail(history$mark, next_history), k]) ) )
    }
    mark_prob = exp(loglik - logSumExp(loglik))
    
    if(sum(mark_prob) == 0) mark_prob = gamma_mat[history$mark[n], ]
    
    next_mark = sample(1:M, 1, prob = mark_prob)
    
    history = rbind(history, list(next_time, next_zone, next_mark, next_history))
    n <- n + 1
  }
  
  return(history[history$time < Tend,])
}

game_simulator_teams <- function(par.sample, time.params, zone.params, meta.g, history, Tend, M){
  
  beta_zone_1 <-  beta_zone_2 <-  beta_zone_3 <- matrix(100, M, M)
  theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
  
  beta_def_raw  = par.sample[1:B1]
  beta_mid_raw  = par.sample[(B1+1):(2*B1)]
  beta_atk_raw  = par.sample[(2*B1+1):(3*B1)]
  theta_def_raw  = par.sample[(3*B1+1):(3*B1+2*B1)]
  theta_mid_raw  = par.sample[(3*B1+2*B1+1):(3*B1+4*B1)]
  theta_atk_raw  = par.sample[(3*B1+4*B1+1):(3*B1+6*B1)]
  
  for (i in 1:B1){
    beta_zone_1[indx_def[i, 1], indx_def[i, 2]] = beta_def_raw[i]
    theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i]
    if(indx_def[i, 1] < (M/2 + 1)){
      beta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }else{
      beta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
      if(indx_def[i, 2] < (M/2)){
        theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
      }
    }
    
    beta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = beta_mid_raw[i]
    theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i]
    if(indx_mid[i, 1] < (M/2 + 1)){
      beta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      } 
    }else{
      beta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
      if(indx_mid[i, 2] < (M/2)){
        theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
      }
    }
    
    beta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = beta_atk_raw[i]
    theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i]
    if(indx_atk[i, 1] < (M/2 + 1)){
      beta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }else{
      beta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
      if(indx_atk[i, 2] < (M/2)){
        theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
      }
    }
  }
  
  beta_list <- list(beta_zone_1, beta_zone_2, beta_zone_3)
  theta <- list(theta_zone_1, theta_zone_2, theta_zone_3)
  
  delta = matrix(par.sample[(9*B1+1):(9*B1+(3*M))], nrow = 3)
  alpha = par.sample[(9*B1+(3*M)+1)]
  
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
  
  shape = time.params[1:15]
  rate = time.params[16:30]
  
  # Initialization
  mark.history <- tail(history, 10)
  mark.n <- nrow(mark.history)
  
  n <- nrow(history)
  last_mark  <- tail(history$mark, 1)
  last_zone  <- tail(history$zone, 1)
  last_score  <- tail(history$score, 1)
  last_history <- tail(history$history, 1)
  
  next_zone <- sample(1:3, 1, prob = 
                        unlist(zone.params[zone.params$p_zone == last_zone & zone.params$p_mark == last_mark, 3:5]))
  
  next_time <- tail(history$time, 1) + rgamma(1, rep(shape, 2)[last_mark], rep(rate, 2)[last_mark])
  
  #next_history = ifelse(last_mark %in% c(10:15, 25:30), 1, min(last_history+1, 5))
  next_history = ifelse(last_mark %in% c(10:15, 25:30), 1, min(last_history+1, 10))
  
  gamma_mat = gamma_list[[next_zone]]
  beta_mat = beta_list[[next_zone]]
  
  loglik <- rep(0, M)
  
  for(k in 1:M){
    loglik[k] = logSumExp(c(log(delta[next_zone, k]),
                            alpha - beta_mat[tail(history$mark, next_history), k] * (next_time - tail(history$time, next_history) ) + 
                              log(gamma_mat[tail(history$mark, next_history), k]) ) )
  }
  mark_prob = exp(loglik - logSumExp(loglik))
  
  if(sum(mark_prob) == 0) mark_prob = gamma_mat[history$mark[n], ]
  
  next_mark = sample(1:M, 1, prob = mark_prob)
  
  history = rbind(history, list(next_time, next_zone, next_mark, next_history))
  n <- n + 1
  
  while(next_time < Tend){
    
    last_mark  <- tail(history$mark, 1)
    last_zone  <- tail(history$zone, 1)
    last_score  <- tail(history$score, 1)
    last_history <- tail(history$history, 1)
    
    next_zone <- sample(1:3, 1, prob = 
                          unlist(zone.params[zone.params$p_zone == last_zone & zone.params$p_mark == last_mark, 3:5]))
    
    next_time <- tail(history$time, 1) + rgamma(1, rep(shape, 2)[last_mark], rep(rate, 2)[last_mark])
    
    #next_history = ifelse(last_mark %in% c(10:15, 25:30), 1, min(last_history+1, 5))
    next_history = ifelse(last_mark %in% c(10:15, 25:30), 1, min(last_history+1, 10))
    
    gamma_mat = gamma_list[[next_zone]]
    beta_mat = beta_list[[next_zone]]
    
    loglik <- rep(0, M)
    
    for(k in 1:M){
      loglik[k] = logSumExp(c(log(delta[next_zone, k]),
                              alpha - beta_mat[tail(history$mark, next_history), k] * (next_time - tail(history$time, next_history) ) + 
                                log(gamma_mat[tail(history$mark, next_history), k]) ) )
    }
    mark_prob = exp(loglik - logSumExp(loglik))
    
    if(sum(mark_prob) == 0) mark_prob = gamma_mat[history$mark[n], ]
    
    next_mark = sample(1:M, 1, prob = mark_prob)
    
    history = rbind(history, list(next_time, next_zone, next_mark, next_history))
    n <- n + 1
  }
  
  return(history[history$time < Tend,])
}

baseline_simulator <- function(base.par, history, Tend, M){
  
  # Initialization
  t <- tail(history$time, 1)
  n <- floor(2*max(base.par)*(Tend - t))
  
  for(i in 1:M){
    for(j in 1:3){
      times <- t + na.omit(cumsum(rexp(n, base.par[i, j])))
      times <- times[times < Tend]
      
      if(length(times) > 0){
        history = rbind(history, cbind(time = times, zone_s = j, mark = i, history = 0))
      }
    }
  }
  history = history[order(history$time),]
  return(history)
}

fbar <- function(u, x, Pmat){
  f <- sapply(1:30, FUN = function(z){
    Pxmin = ecdf(Pmat[z,])(x[z]-1)
    Px = ecdf(Pmat[z,])(x[z])
    ret <- (u - Pxmin)/(Px - Pxmin)
    ret <- ifelse(ret < 0, 0, ifelse(ret > 1, 1, ret))
  })
  mean(f)
}

# Markov location model
model.likelihood.zones <- function(params, period.data, M = 30, cl){
  
  period.data[, p_mark := shift(mark, 1, type = 'lag'), by = 'id']
  period.data[, p_zone := shift(zone_s, 1, type = 'lag'), by = 'id']
  
  liks = apply(params, 1, FUN = function(z) {
    vec <- rep(0, nrow(period.data) - 1)
    for(i in 2:nrow(period.data)){
      prob <- z[(1 + 30*(period.data$p_zone[i]-1)):(30*period.data$p_zone[i]),]
      vec[i-1] <- prob[period.data$p_mark[i], period.data$zone_s[i]]
    }
    vec
  })
  
  loglik <- log(liks)
  
  return(loglik)
}

# Gamma times model
model.likelihood.times <- function(params, period.data, M = 30, cl){
  
  clusterExport(cl, list("period.data", "params", "M"), envir=environment())
  
  period.likl.array <- parApply(cl, params, 1, FUN = function(par.sample){
    
    shape = rep(par.sample[(M*(M-1)+M+3):(M*(M-1)+M+17)], 2)
    rate  = rep(par.sample[(M*(M-1)+M+18):(M*(M-1)+M+32)], 2)
    
    n <- nrow(period.data)
    time_vec <- period.data$time
    mark_vec <- period.data$mark
    
    ll <- dgamma(diff(time_vec), shape = shape[head(mark_vec, -1)], rate = rate[head(mark_vec, -1)], log = T)
  })
}

# Simple beta model
model.likelihood.marks.simple <- function(params, period.data, M = 30, cl){
  
  clusterExport(cl, list("period.data", "params", "M"), envir=environment())
  
  period.likl.array <- parApply(cl, params, 1, FUN = function(par.sample){
    
    delta = par.sample[1:M]
    alpha = par.sample[M+1]
    beta  = par.sample[M+2]
    theta = matrix(par.sample[(M+3):(M+2+M*(M-1))], nrow = M, ncol = M-1)
    sigma_gamma = par.sample[M*M+M+3]
    theta = theta*sigma_gamma
    gamma <- exp(theta)/(1+rowSums(exp(theta)))
    gamma <- cbind(gamma, 1 - rowSums(gamma))
    
    n <- nrow(period.data)
    time_vec <- period.data$time
    mark_vec <- period.data$mark
    hist_vec <- period.data$history
    ll <- rep(0, n-1)
    
    for(i in 2:n){
      
      ll[i-1] = 
        logSumExp(c(log(delta[mark_vec[i]]),
                    alpha - beta * (time_vec[i] - time_vec[(i-hist_vec[i]):(i-1)] ) + 
                      log(gamma[mark_vec[(i-hist_vec[i]):(i-1)], mark_vec[i]]) ) ) - 
        logSumExp(c(0,
                    alpha - beta * (time_vec[i] - time_vec[(i-hist_vec[i]):(i-1)] ) ))
    }
    return(ll)
  })
}

# Vector beta model
model.likelihood.marks.vector <- function(params, period.data, M = 30, cl){
  
  clusterExport(cl, list("period.data", "params", "M"), envir=environment())
  
  period.likl.array <- parApply(cl, params, 1, FUN = function(par.sample){
    
    delta = par.sample[1:M]
    alpha = par.sample[M+1]
    beta  = rep(par.sample[(M+2):(M+16)],2)
    theta = matrix(par.sample[(M+17):(M+16+M*(M-1))], nrow = M, ncol = M-1)
    sigma_gamma = par.sample[M*M+M+17]
    theta = theta*sigma_gamma
    gamma <- exp(theta)/(1+rowSums(exp(theta)))
    gamma <- cbind(gamma, 1 - rowSums(gamma))
    
    n <- nrow(period.data)
    time_vec <- period.data$time
    mark_vec <- period.data$mark
    hist_vec <- period.data$history
    ll <- rep(0, n-1)
    
    for(i in 2:n){
      
      ll[i-1] = 
        logSumExp(c(log(delta[mark_vec[i]]),
                    alpha - beta[mark_vec[(i-hist_vec[i]):(i-1)]] * (time_vec[i] - time_vec[(i-hist_vec[i]):(i-1)] ) + 
                      log(gamma[mark_vec[(i-hist_vec[i]):(i-1)], mark_vec[i]]) ) ) - 
        logSumExp(c(0,
                    alpha - beta[mark_vec[(i-hist_vec[i]):(i-1)]] * (time_vec[i] - time_vec[(i-hist_vec[i]):(i-1)] ) ))
    }
    return(ll)
  })
}

# Matrix beta model
model.likelihood.marks <- function(params, period.data, M = 30, cl){
  
  clusterExport(cl, list("period.data", "params", "M", "indx_def", "indx_mid", "indx_atk", "B1"), envir=environment())
  
  period.likl.array <- parApply(cl, params, 1, FUN = function(par.sample){
    
    beta_zone_1 <-  beta_zone_2 <-  beta_zone_3 <- matrix(100, M, M)
    theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
    
    beta_def_raw  = par.sample[1:B1]
    beta_mid_raw  = par.sample[(B1+1):(2*B1)]
    beta_atk_raw  = par.sample[(2*B1+1):(3*B1)]
    theta_def_raw  = par.sample[(3*B1+1):(5*B1)]
    theta_mid_raw  = par.sample[(5*B1+1):(7*B1)]
    theta_atk_raw  = par.sample[(7*B1+1):(9*B1)]
    
    for (i in 1:B1){
      beta_zone_1[indx_def[i, 1], indx_def[i, 2]] = beta_def_raw[i]
      theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i]
      if(indx_def[i, 1] < (M/2 + 1)){
        beta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
        if(indx_def[i, 2] < (M/2)){
          theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
        }
      }else{
        beta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
        if(indx_def[i, 2] < (M/2)){
          theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
        }
      }
      
      beta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = beta_mid_raw[i]
      theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i]
      if(indx_mid[i, 1] < (M/2 + 1)){
        beta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
        if(indx_mid[i, 2] < (M/2)){
          theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
        } 
      }else{
        beta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
        if(indx_mid[i, 2] < (M/2)){
          theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
        }
      }
      
      beta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = beta_atk_raw[i]
      theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i]
      if(indx_atk[i, 1] < (M/2 + 1)){
        beta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
        if(indx_atk[i, 2] < (M/2)){
          theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
        }
      }else{
        beta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
        if(indx_atk[i, 2] < (M/2)){
          theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
        }
      }
    }
    
    beta_list <- list(beta_zone_1, beta_zone_2, beta_zone_3)
    theta <- list(theta_zone_1, theta_zone_2, theta_zone_3)
    
    delta = matrix(par.sample[(9*B1+1):(9*B1+(3*M))], nrow = 3)
    alpha = par.sample[(9*B1+(3*M)+1)]
    
    gamma_list <- rep(list(matrix(0, nrow = M, ncol = M-1)),3)
    
    for(i in 1:3){
      gamma_iter <- theta[[i]]
      
      gamma_iter <- exp(gamma_iter)/(1+rowSums(exp(gamma_iter)))
      gamma_iter <- cbind(gamma_iter, 1 - rowSums(gamma_iter))
      
      gamma_list[[i]] <- gamma_iter
    }
    
    n <- nrow(period.data)
    time_vec <- period.data$time
    mark_vec <- period.data$mark
    loc_vec <- period.data$zone_s
    hist_vec <- period.data$history
    loglik <- rep(0, M)
    ll <- rep(0, n-1)
    
    for(i in 2:n){
      beta_mat = beta_list[[loc_vec[i]]]
      gamma_mat = gamma_list[[loc_vec[i]]]
      
      for(k in 1:M){
        loglik[k] = logSumExp(c(log(delta[loc_vec[i], k]),
                                alpha - beta_mat[mark_vec[(i-hist_vec[i]):(i-1)], k] * (time_vec[i] - time_vec[(i-hist_vec[i]):(i-1)] ) + 
                                  log(gamma_mat[mark_vec[(i-hist_vec[i]):(i-1)], k]) ) )
      }
      ll[i-1] = loglik[mark_vec[i]] - logSumExp(loglik)
    }
    return(ll)
  })
}

# Matrix beta with teams
model.likelihood.marks.teams <- function(params, period.data, meta.g, M = 30, cl){
  
  clusterExport(cl, list("period.data", "params", "M", "indx_def", "indx_mid", "indx_atk", "B1", "meta.g"), envir=environment())
  
  period.likl.array <- parApply(cl, params, 1, FUN = function(par.sample){
    
    beta_zone_1 <-  beta_zone_2 <-  beta_zone_3 <- matrix(100, M, M)
    theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
    
    beta_def_raw  = par.sample[1:B1]
    beta_mid_raw  = par.sample[(B1+1):(2*B1)]
    beta_atk_raw  = par.sample[(2*B1+1):(3*B1)]
    theta_def_raw  = par.sample[(3*B1+1):(3*B1+2*B1)]
    theta_mid_raw  = par.sample[(3*B1+2*B1+1):(3*B1+4*B1)]
    theta_atk_raw  = par.sample[(3*B1+4*B1+1):(3*B1+6*B1)]
    
    for (i in 1:B1){
      beta_zone_1[indx_def[i, 1], indx_def[i, 2]] = beta_def_raw[i]
      theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i]
      if(indx_def[i, 1] < (M/2 + 1)){
        beta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
        if(indx_def[i, 2] < (M/2)){
          theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
        }
      }else{
        beta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
        if(indx_def[i, 2] < (M/2)){
          theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
        }
      }
      
      beta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = beta_mid_raw[i]
      theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i]
      if(indx_mid[i, 1] < (M/2 + 1)){
        beta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
        if(indx_mid[i, 2] < (M/2)){
          theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
        } 
      }else{
        beta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
        if(indx_mid[i, 2] < (M/2)){
          theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
        }
      }
      
      beta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = beta_atk_raw[i]
      theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i]
      if(indx_atk[i, 1] < (M/2 + 1)){
        beta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
        if(indx_atk[i, 2] < (M/2)){
          theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
        }
      }else{
        beta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
        if(indx_atk[i, 2] < (M/2)){
          theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
        }
      }
    }
    
    beta_list <- list(beta_zone_1, beta_zone_2, beta_zone_3)
    theta <- list(theta_zone_1, theta_zone_2, theta_zone_3)
    
    delta = matrix(par.sample[(9*B1+1):(9*B1+(3*M))], nrow = 3)
    alpha = par.sample[(9*B1+(3*M)+1)]
    
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
    
    n <- nrow(period.data)
    time_vec <- period.data$time
    mark_vec <- period.data$mark
    loc_vec <- period.data$zone_s
    hist_vec <- period.data$history
    loglik <- rep(0, M)
    ll <- rep(0, n-1)
    
    for(i in 2:n){
      beta_mat = beta_list[[loc_vec[i]]]
      gamma_mat = gamma_list[[loc_vec[i]]]
      
      for(k in 1:M){
        loglik[k] = logSumExp(c(log(delta[loc_vec[i], k]),
                                alpha - beta_mat[mark_vec[(i-hist_vec[i]):(i-1)], k] * (time_vec[i] - time_vec[(i-hist_vec[i]):(i-1)] ) + 
                                  log(gamma_mat[mark_vec[(i-hist_vec[i]):(i-1)], k]) ) )
      }
      ll[i-1] = loglik[mark_vec[i]] - logSumExp(loglik)
    }
    return(ll)
  })
}

# Matrix beta with teams delta
model.likelihood.marks.teams.delta <- function(params, period.data, meta.g, M = 30, cl){
  
  clusterExport(cl, list("period.data", "params", "M", "indx_def", "indx_mid", "indx_atk", "B1", "meta.g"), envir=environment())
  
  period.likl.array <- parApply(cl, params, 1, FUN = function(par.sample){
    
    beta_zone_1 <-  beta_zone_2 <-  beta_zone_3 <- matrix(100, M, M)
    theta_zone_1 <- theta_zone_2 <- theta_zone_3 <- matrix(-100, M, M-1)
    
    beta_def_raw  = par.sample[1:B1]
    beta_mid_raw  = par.sample[(B1+1):(2*B1)]
    beta_atk_raw  = par.sample[(2*B1+1):(3*B1)]
    theta_def_raw  = par.sample[(3*B1+1):(3*B1+2*B1)]
    theta_mid_raw  = par.sample[(3*B1+2*B1+1):(3*B1+4*B1)]
    theta_atk_raw  = par.sample[(3*B1+4*B1+1):(3*B1+6*B1)]
    
    for (i in 1:B1){
      beta_zone_1[indx_def[i, 1], indx_def[i, 2]] = beta_def_raw[i]
      theta_zone_1[indx_def[i, 1], indx_def[i, 2]] = theta_def_raw[i]
      if(indx_def[i, 1] < (M/2 + 1)){
        beta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
        if(indx_def[i, 2] < (M/2)){
          theta_zone_3[indx_def[i, 1]+(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
        }
      }else{
        beta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = beta_def_raw[i]
        if(indx_def[i, 2] < (M/2)){
          theta_zone_3[indx_def[i, 1]-(M/2), indx_def[i, 2]+(M/2)] = theta_def_raw[B1+i]
        }
      }
      
      beta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = beta_mid_raw[i]
      theta_zone_2[indx_mid[i, 1], indx_mid[i, 2]] = theta_mid_raw[i]
      if(indx_mid[i, 1] < (M/2 + 1)){
        beta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
        if(indx_mid[i, 2] < (M/2)){
          theta_zone_2[indx_mid[i, 1]+(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
        } 
      }else{
        beta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = beta_mid_raw[i]
        if(indx_mid[i, 2] < (M/2)){
          theta_zone_2[indx_mid[i, 1]-(M/2), indx_mid[i, 2]+(M/2)] = theta_mid_raw[B1+i]
        }
      }
      
      beta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = beta_atk_raw[i]
      theta_zone_3[indx_atk[i, 1], indx_atk[i, 2]] = theta_atk_raw[i]
      if(indx_atk[i, 1] < (M/2 + 1)){
        beta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
        if(indx_atk[i, 2] < (M/2)){
          theta_zone_1[indx_atk[i, 1]+(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
        }
      }else{
        beta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = beta_atk_raw[i]
        if(indx_atk[i, 2] < (M/2)){
          theta_zone_1[indx_atk[i, 1]-(M/2), indx_atk[i, 2]+(M/2)] = theta_atk_raw[B1+i]
        }
      }
    }
    
    beta_list <- list(beta_zone_1, beta_zone_2, beta_zone_3)
    theta <- list(theta_zone_1, theta_zone_2, theta_zone_3)
    
    delta1 = matrix(par.sample[(9*B1+1):(9*B1+M)], nrow = 1)
    delta2 = matrix(par.sample[(9*B1+1+M+14+19*(M-1)+1):(9*B1+1+M+14+19*(M-1)+14)], nrow = 1)
    delta2 = c(delta2, 0.5 - sum(delta2))
    delta  = rbind(delta1, c(delta2, delta2), c(delta1[16:30], delta1[1:15]))
    alpha  = par.sample[(9*B1+M+15)]
    
    mu_atk  <- rbind(matrix(par.sample[(9*B1+1+M+14+1):(9*B1+1+M+14+19*(M-1))], nrow = 19), rep(0, M-1))
    
    mu.game <- c(mu_atk[meta.g[2], 1:(M/2)], mu_atk[meta.g[3], -c(1:(M/2))])
    
    gamma =  matrix(rep(mu.game, each = M), nrow = M, ncol = M-1)
    gamma_list <- rep(list(gamma),3)
    
    for(i in 1:3){
      
      gamma_iter <- theta[[i]] + gamma
      
      gamma_iter <- exp(gamma_iter)/(1+rowSums(exp(gamma_iter)))
      gamma_iter <- cbind(gamma_iter, 1 - rowSums(gamma_iter))
      
      gamma_list[[i]] <- gamma_iter
    }
    
    n <- nrow(period.data)
    time_vec <- period.data$time
    mark_vec <- period.data$mark
    loc_vec <- period.data$zone_s
    hist_vec <- period.data$history
    loglik <- rep(0, M)
    ll <- rep(0, n-1)
    
    for(i in 2:n){
      beta_mat = beta_list[[loc_vec[i]]]
      gamma_mat = gamma_list[[loc_vec[i]]]
      
      for(k in 1:M){
        loglik[k] = logSumExp(c(log(delta[loc_vec[i], k]),
                                alpha - beta_mat[mark_vec[(i-hist_vec[i]):(i-1)], k] * (time_vec[i] - time_vec[(i-hist_vec[i]):(i-1)] ) + 
                                  log(gamma_mat[mark_vec[(i-hist_vec[i]):(i-1)], k]) ) )
      }
      ll[i-1] = loglik[mark_vec[i]] - logSumExp(loglik)
    }
    return(ll)
  })
}

# Homo Baseline model
base.likelihood<- function(params, period.data, M = 30){
  
  counts = as.numeric(table(factor(period.data$mark, levels = 1:M), period.data$zone_s))
  totaltime = max(period.data$time)
  
  loglik <- apply(params, 1, FUN = function(z) sum(counts*log(z) - totaltime*z))
  
  return(loglik)
}

# Smart baseline model
smartbase.likelihood.marks <- function(params, base.params, period.data, M = 30){
  
  period.data[, p_mark := shift(mark, 1, type = 'lag')]
  
  liks = apply(params, 1, FUN = function(z) {
    z[cbind(period.data$p_mark[-1], period.data$mark[-1])]
  })
  
  liks[which(rowSums(liks) == 0), ] = base.params[period.data$mark[which(rowSums(liks) == 0)+1]]
  loglik <- log(liks)
  
  return(loglik)
}

# ROC plot
ggrocs <- function(rocs, breaks = seq(0,1,0.1), legendTitel = "Legend (AUC)") {
    require(plyr)
    RocVals <- plyr::ldply(names(rocs), function(rocName) {
      
      data.frame(cbind(coords(rocs[[rocName]], x = seq(0, 1, 0.01), input="threshold", ret = c("threshold", "tpr", "fpr"), transpose = FALSE),
            names = rep(rocName, 101)))
    })
    RocVals$names <- factor(RocVals$names, levels = c('Model', 'MA_5', 'MA_10', 'MA_15'))
    RocVals <- RocVals[order(RocVals$names, RocVals$fpr, RocVals$tpr),]
    aucs <- sapply(rocs, "[[", "auc")
    
    rocPlot <- ggplot(RocVals, aes(x = fpr, y = tpr, colour = names)) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), alpha = 0.5, colour = "gray") + 
      geom_point() +
      geom_line() +
      scale_x_continuous(name = "False Positive Rate (1 - Specificity)", limits = c(0,1), breaks = breaks) + 
      scale_y_continuous(name = "True Positive Rate (Sensitivity)", limits = c(0,1), breaks = breaks) +
      theme_bw() + 
      coord_equal() + 
      scale_color_manual(labels = paste0(names(rocs), ' (', round(aucs,2) , ')'), breaks = unique(RocVals$names), 
                         values = c('red', 'blue', 'green', 'orange')) +
      guides(colour = guide_legend(legendTitel)) +
      theme(axis.ticks = element_line(color = "grey80"))
    rocPlot
}
