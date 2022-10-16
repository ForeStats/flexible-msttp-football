## Distributed as part of the supporting materials for the manuscript
## "Flexible marked spatio-temporal point processes with applications to event sequences from association football"
##
## Author: Santhosh Narayanan
## Date: 16 Oct 2022
## Licence: GPL 2 or greater
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! Provided "as is"
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

## Helper functions for Section 3

simulate_multi_hawkes <- function(params, n, d) {
  #Parameters
  A       <- params[[1]]
  beta    <- params[[2]]
  mu      <- params[[3]]
  
  #Initialization
  U       <- runif(d)
  t       <- min(-log(U) / mu)
  l       <- which(-log(U) / mu == t)
  events  <- data.frame(time = t, eid = l)
  lambdak <- mu + A[l, ] * beta
  k       <- 1
  
  while (k < n) {
    U1       <- runif(d)
    D        <- 1 + (beta * log(U1) / (lambdak - mu))
    D[D < 0] <- 0
    S1       <- -log(D) / beta
    
    U2 <- runif(d)
    S2 <- -log(U2) / mu
    
    S  <- colMins(rbind(S1, S2))
    
    W  <- min(S)
    l  <- which(S == W)
    events  <- rbind(events, c(events$time[k] + W, l))
    lambdak <- mu + A[l, ] * beta + (lambdak - mu) * exp(-beta * W)
    k  <- k + 1
  }
  
  return(events)
}

# Log-likelihood function
loglik <- function(params, events, endT, d) {
  A       <- matrix(params[1:(d * d)], nrow = d)
  B       <- params[(d * d + 1):(d * (d + 1))]
  mu      <- params[(d * (d + 1) + 1):(d * (d + 2))]
  
  n       <- nrow(events)
  
  term_1 <- -sum(exp(mu) * endT)
  
  term_2 <- 0
  
  if (n > 0) {
    tdiff <- endT - events$time
    Au <- A[events$eid,]
    term_2 <-
      sum(exp(Au) * t(exp(-exp(B) %*% t(tdiff)) - 1), na.rm = T)
  }
  
  K <- 0
  if (n > 1) {
    K <- c(0, sapply(2:n, function(z) {
      Bz <- B[events$eid[z]]
      Az <- A[events$eid[1:(z - 1)], events$eid[z]]
      tdiff <- events$time[z] - events$time[1:(z - 1)]
      sum(exp(Az) * exp(Bz) * exp(-exp(Bz) * (tdiff)), na.rm = T)
    }))
  }
  lambdak <- mu[events$eid]
  
  term_3 <- sum(log(exp(lambdak) + K))
  
  loglik <- -term_1 - term_2 - term_3
  
  return(loglik)
}

# gradient
gradloglik <- function(params, events, endT, d) {
  A       <- matrix(params[1:(d * d)], nrow = d)
  B       <- params[(d * d + 1):(d * (d + 1))]
  mu      <- params[(d * (d + 1) + 1):(d * (d + 2))]
  
  n       <- nrow(events)
  
  # Mu
  mu_term <- rep(0, d)
  for (u in 1:d) {
    mu_d <- mu[u]
    term1 <- -exp(mu_d) * endT
    
    term2 <- 0
    subset <- events[events$eid == u, ]
    nu <- nrow(subset)
    
    if (nu > 0) {
      K <- sapply(1:nu, function(z) {
        t <- subset$time[z]
        Bz <- B[subset$eid[z]]
        Az <- A[events$eid[events$time < t], subset$eid[z]]
        tdiff <- t - events$time[events$time < t]
        sum(exp(Az) * exp(Bz) * exp(-exp(Bz) * (tdiff)), na.rm = T)
      })
      
      term2 <- sum(exp(mu_d) / (exp(mu_d) + K))
    }
    mu_term[u] <- term1 + term2
    
  }
  
  # beta
  beta_term <- rep(0, d)
  for (u in 1:d) {
    Au <- A[events$eid, u]
    Bu <- B[u]
    tdiff <- endT - events$time
    term_1 <-
      sum(exp(Au) * (tdiff) * exp(-exp(Bu) * tdiff) * exp(Bu), na.rm = T)
    
    term_2 <- 0
    subset <- events[events$eid == u, ]
    nu <- nrow(subset)
    if (nu > 0) {
      K1 <- sapply(1:nu, function(z) {
        t <- subset$time[z]
        Bz <- B[subset$eid[z]]
        Az <- A[events$eid[events$time < t], subset$eid[z]]
        tdiff <- t - events$time[events$time < t]
        sum(exp(Az) * exp(Bz) * exp(-exp(Bz) * tdiff), na.rm = T)
      })
      K2 <- sapply(1:nu, function(z) {
        t <- subset$time[z]
        Bz <- B[subset$eid[z]]
        Az <- A[events$eid[events$time < t], subset$eid[z]]
        tdiff <- t - events$time[events$time < t]
        sum(exp(Az) * exp(Bz) * exp(-exp(Bz) * tdiff) * (1 - exp(Bz) * tdiff),
            na.rm = T)
      })
      mu_d <- mu[u]
      term_2 <- sum(K2 / (exp(mu_d) + K1), na.rm = T)
    }
    
    beta_term[u] <- -term_1 + term_2
  }
  
  # Alpha
  alpha_term <- matrix(0, nrow = d, ncol = d)
  for (i in 1:d) {
    for (j in 1:d) {
      Aij <- A[i, j]
      Bj <- B[j]
      
      subset <- events[events$eid == i, ]
      
      term_1 <- 0
      if (nrow(subset) > 0) {
        tdiff <- endT - subset$time
        term_1 <-
          sum(exp(Aij) * (exp(-exp(Bj) * (tdiff)) - 1), na.rm = T)
      }
      
      subset1 <- events[events$eid == j, ]
      nu <- nrow(subset1)
      
      term_2 <- 0
      
      if (nu > 0 & nrow(subset) > 0) {
        K1 <- sapply(1:nu, function(z) {
          t <- subset1$time[z]
          Az <- A[events$eid[events$time < t], j]
          tdiff <- t - events$time[events$time < t]
          sum(exp(Az) * exp(Bj) * exp(-exp(Bj) * (tdiff)), na.rm = T)
        })
        
        K2 <- sapply(1:nu, function(z) {
          t <- subset1$time[z]
          tdiff <- t - subset$time[subset$time < t]
          sum(exp(Aij) * exp(Bj) * exp(-exp(Bj) * tdiff), na.rm = T)
        })
        
        mu_d <- mu[j]
        
        term_2 <- sum(K2 / (exp(mu_d) + K1), na.rm = T)
      }
      
      alpha_term[i, j] <- term_1 +  term_2
    }
  }
  
  return(-c(alpha_term, beta_term, mu_term))
}

# Cumulated Log-likelihood function
cumloglik <- function(params, events, endtimes, n, d) {
  cat(".")
  
  cumloglik <- foreach(sample = 1:n,
                       .combine = sum,
                       .multicombine = TRUE)  %dopar%  {
                         loglik(params, events = events[which(events$Segment == sample),
                                                        c("time", "eid")], endT = endtimes[sample], d)
                       }
  
  grad <- foreach(sample = 1:n,
                  .combine = rbind,
                  .multicombine = TRUE)  %dopar%  {
                    gradloglik(params, events = events[which(events$Segment == sample),
                                                       c("time", "eid")], endT = endtimes[sample], d)
                  }
  
  grad <-  colSums(grad, na.rm = T)
  attr(cumloglik, "gradient") <- grad
  
  return(cumloglik)
}

InHomoK <- function(s, events, lims) {
  N <- length(events)
  
  outersum <- sapply(
    events,
    FUN = function(z) {
      xdash <- events[!events %in% z]
      
      Ks <-
        ifelse(abs(xdash - z) <= s, ifelse(abs(xdash - z) <= min(z, lims[2] - z), 1, 2), 0)
      
      sum(Ks)
    }
  )
  
  sum(outersum) * (lims[2] - lims[1]) / (N ^ 2)
}

simulate_uni_hawkes <- function(params, n) {
  #Parameters
  mu    <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  #Initialization
  U       <- runif(1)
  t       <- min(-log(U) / mu)
  events  <- t
  lambdak <- mu + alpha * beta
  k       <- 1
  
  while (k < n) {
    D  <- 1 + (beta * log(runif(1)) / (lambdak - mu))
    D[D < 0] <- 0
    S1 <- -log(D) / beta
    
    S2 <- -log(runif(1)) / mu
    
    W <- min(S1, S2)
    events  <- c(events, events[k] + W)
    lambdak <- mu + alpha * beta + (lambdak - mu) * exp(-beta * W)
    
    k = k + 1
  }
  
  return(events)
}

# Log-likelihood function
uniloglik <- function(params, events, endT) {
  mu    <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  n       <- length(events)
  
  term_1 <- -sum(exp(mu) * endT)
  
  term_2 <- 0
  
  if (n > 0) {
    tdiff <- endT - events
    term_2 <-
      sum(exp(alpha) * t(exp(-exp(beta) %*% t(tdiff)) - 1), na.rm = T)
  }
  
  K <- 0
  if (n > 1) {
    K <- c(0, sapply(2:n, function(z) {
      tdiff <- events[z] - events[1:(z - 1)]
      sum(exp(alpha) * exp(beta) * exp(-exp(beta) * (tdiff)), na.rm = T)
    }))
  }
  
  term_3 <- sum(log(exp(mu) + K))
  
  loglik <- -term_1 - term_2 - term_3
  
  return(loglik)
}

# Fixed Log-likelihood function
fixedloglik <- function(params, events, endT, alpha, beta) {
  mu     <- params[1]
  
  n      <- length(events)
  
  term_1 <- -sum(exp(mu) * endT)
  
  term_2 <- 0
  
  if (n > 0) {
    tdiff  <- endT - events
    term_2 <-
      sum(exp(alpha) * t(exp(-exp(beta) %*% t(tdiff)) - 1), na.rm = T)
  }
  
  K <- 0
  if (n > 1) {
    K <- c(0, sapply(2:n, function(z) {
      tdiff <- events[z] - events[1:(z - 1)]
      sum(exp(alpha) * exp(beta) * exp(-exp(beta) * (tdiff)), na.rm = T)
    }))
  }
  
  term_3 <- sum(log(exp(mu) + K))
  
  loglik <- -term_1 - term_2 - term_3
  
  return(loglik)
}
