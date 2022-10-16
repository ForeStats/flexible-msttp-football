## Distributed as part of the supporting materials for the manuscript
## "Flexible marked spatio-temporal point processes with applications to event sequences from association football"
##
## Author: Santhosh Narayanan
## Date: 16 Oct 2022
## Licence: GPL 2 or greater
## NOT A POLISHED PIECE OF PUBLIC-USE SOFTWARE! Provided "as is"
## NO WARRANTY OF FITNESS FOR ANY PURPOSE!

## This script generates Figure 4,5 in Section 3

library(data.table)
library(doParallel)
library(matrixStats)
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)

#Set seed
set.seed(001)

# Load the data after Data Wrangling Script
# load("matchData.RData")
source("helpers.R")

# Set sample size for experiment
n.events <- 10000

game.data  <-
  data_season_1[game %in% 1:8, c('game', 'period', 'time')]
rm(
  list = c(
    "data_season_1",
    "data_season_2",
    "df_season_1",
    "df_season_2",
    "mark_levels"
  )
)
game_times <- diff(game.data$time)
game_times <- game_times + rnorm(length(game_times), 0, 0.4)
game_times <- game_times[game_times > 0]
game_times <- game_times[1:n.events]

# Fit Hawkes Model
hawkes.data <- copy(game.data)
hawkes.data <-
  hawkes.data[, Segment := sequence(.N), by = c('game', 'period')]
hawkes.data[Segment > 1, Segment := 0]
hawkes.data[, Segment := cumsum(Segment)]
hawkes.data[, game := NULL]
hawkes.data[, period := NULL]
hawkes.data[, time := time - min(time), by = 'Segment']
hawkes.data[, tdiff := c(0, diff(time)), by = 'Segment']
hawkes.data[tdiff == 0, tdiff := 0.1]
hawkes.data[, time := cumsum(tdiff), by = 'Segment']
hawkes.data[, eid := 1]
hawkes.data[, tdiff := NULL]

endT <- unlist(hawkes.data[, tail(time, 1), by = "Segment"][, 2])

cl   <- makePSOCKcluster(detectCores() - 1)
registerDoParallel(cl)
clusterExport(cl, varlist = c('loglik', 'gradloglik'))

par <-  log(c(0.0001, 0.0001, 0.382))
estimates <-
  nlm(
    p = par,
    f = cumloglik,
    events = data.frame(hawkes.data[, 1:3]),
    endtimes = endT,
    n = 16,
    d = 1,
    gradtol = 1e-5,
    print.level = 2
  )

stopCluster(cl)

par.est       <- exp(estimates$estimate)
params0       <- list(matrix(par.est[1]), par.est[2], par.est[3])
hawkes.events <- data.table()
for (i in 1:16) {
  hawkes.events <-
    rbind(hawkes.events, cbind(Segment = i,
                               simulate_multi_hawkes(params0, n = 626, d = 1)))
}
hawkes_times <- diff(hawkes.events$time)
hawkes_times <- hawkes_times[hawkes_times > 0]
  
# Poisson Fit
pois_times  <- rexp(n.events, rate = 1 / mean(game_times))

# Gamma fit
rate        <- mean(game_times) / var(game_times)
shape       <- (mean(game_times) ^ 2) / var(game_times)
gamma_times <- rgamma(n.events, shape = shape, rate = rate)

# K-S test
ks.test(hawkes_times, pois_times)
ks.test(game_times, pois_times)
ks.test(game_times, hawkes_times)
ks.test(game_times, gamma_times)

n    <- n.events
vals <- data.frame(
  Poisson = pois_times[1:n],
  Hawkes  = hawkes_times[1:n],
  Gamma   = gamma_times[1:n],
  Football = game_times[1:n]
)

myColours <- brewer.pal(4,"Dark2")

my.settings <- standard.theme("pdf", color=TRUE)
my.settings$superpose.line$lwd <- 2
my.settings$superpose.line$col <- myColours[1:4]

pdf("figures/ecdf.pdf", width=10, height=8)
ecdfplot(
  ~ Poisson + Hawkes + Gamma + Football,
  data = vals,
  xlab = textGrob('inter arrival time', gp = gpar(
    fontfamily = "CMU Serif", fontsize = 20
  )),
  ylab = textGrob(
    'cumulative distribution function',
    rot = 90,
    gp = gpar(fontfamily = "CMU Serif", fontsize = 20)
  ),
  xlim = c(-1, 10),
  lwd = 1.2,
  auto.key = list(
    text = c(
      'Poisson process',
      'Hawkes process',
      'Gamma process',
      'Observed times'
    ),
    space = "top",
    columns = 2,
    cex = 1.5,
    fontfamily = "CMU Serif"
  ),
  par.settings = my.settings,
  scales = list(tck = c(1, 0)),
  family = "CMU Serif"
) 
dev.off()
embed_fonts("figures/ecdf.pdf", outfile="figures/ecdf.pdf")

## Second order analysis
game.data = game.data[game == 1, c('game', 'period', 'time')]

game_times = diff(game.data$time)
rm(list = ls()[!ls() %in% c('game_times')])
source("helpers.R")

game_times <- game_times + rnorm(length(game_times), 0, 0.2)
game_times <- game_times[game_times > 0]
game_times <- cumsum(game_times)
lims       <- c(0, round(max(game_times) + 0.1, 0))
game_intensity <-  length(game_times) / (lims[2] - lims[1])

game_values <-
  sapply(
    seq(1, 100, 1),
    FUN = function(k)
      InHomoK(s = k, events = game_times, lims = lims) - 2 * k
  )
plot(game_values)

# MLE
par <-  log(c(0.5, 0.5, 0.1))
estimates <-
  nlm(
    p = par,
    f = uniloglik,
    events = game_times,
    endT = lims[2],
    gradtol = 1e-6,
    print.level = 2
  )
par.est.mle <- exp(estimates$estimate)

# fixed MLE
param_df <-
  cbind(rep(par.est.mle[1], 3),
        c(par.est.mle[2], 0.4, 0.8),
        c(par.est.mle[3], 0.01, 0.01))

for(i in 2:3) {
  par <-  log(par.est.mle[1])
  estimates <-
    nlm(
      p = par,
      f = fixedloglik,
      events = game_times,
      endT = lims[2],
      alpha = log(param_df[i, 2]),
      beta = log(param_df[i, 3]),
      gradtol = 1e-6,
      print.level = 2
    )
  param_df[i, 1] <- exp(estimates$estimate)
}

# Simulations
Nsim <- 100
hawkes_arr <- array(NA, dim = c(3, Nsim, 100)) 

for(j in 1:3) {
  for (i in 1:Nsim) {
    hawkes_times <-
      simulate_uni_hawkes(params = unlist(param_df[j, ]),
                          n = 2 * length(game_times))
    hawkes_times <- hawkes_times[hawkes_times < lims[2]]
    
    hawkes_arr[j, i,] <-
      sapply(
        seq(1, 100, 1),
        FUN = function(k)
          InHomoK(
            s = k,
            events = hawkes_times,
            lims = lims
          ) - 2 * k
      )
  }
}

poisson_df <- data.frame()
for (i in 1:Nsim) {
  poisson_times <-
    cumsum(rexp(n = 2 * length(game_times), r = game_intensity))
  poisson_times <- poisson_times[poisson_times < lims[2]]
  
  poisson_values <-
    sapply(
      seq(1, 100, 1),
      FUN = function(k)
        InHomoK(s = k, events = poisson_times, lims = lims) - 2 * k
    )
  poisson_df <- rbind(poisson_df, poisson_values)
}

myColours <- brewer.pal(4, "Dark2")

plot_data <-
  data.frame(x = seq(1, 100, 1), Hawkes = t(hawkes_arr[1, ,]))
plot_data <- melt(plot_data, id = 1)
plot_data$variable <- as.character(plot_data$variable)
plot_data$variable <- 'Hawkes'

#  Boxplot
p1 <- ggplot() +
  geom_point(aes(x = factor(1:100), y = game_values, fill = 'Observed'))  +
  scale_fill_manual(name = NULL,
                    labels = c('Observed times'),
                    values = 'black') +
  geom_boxplot(
    data = plot_data,
    aes(x = factor(x), y = value, colour = 'Hawkes'),
    outlier.shape = NA
  ) +
  scale_colour_manual(
    name = NULL,
    labels = c('Hawkes I'),
    values = c("Hawkes" = myColours[1])
  ) +
  scale_x_discrete(breaks = c(0, 25, 50, 75)) +
  scale_y_continuous(
    breaks = c(-2, 0, 2),
    labels = c("-  2", "0", "2"),
    limits = c(-3, 3)
  ) +
  labs(colour = NULL, fill = NULL) +
  xlab("t") + ylab("K(t) -  2t") +
  theme_bw() +
  theme(
    legend.position = "top",
    text = element_text(
      family = "CMU Serif",
      size = 14,
      vjust = 0.9
    )
  )

plot_data <-
  data.frame(x = seq(1, 100, 1), Hawkes = t(hawkes_arr[2, , ]))
plot_data <- melt(plot_data, id = 1)
plot_data$variable <- as.character(plot_data$variable)
plot_data$variable <- 'Hawkes'

#  Boxplot
p2 <- ggplot() +
  geom_point(aes(x = factor(1:100), y = game_values, fill = 'Observed'))  +
  scale_fill_manual(name = NULL,
                    labels = c('Observed times'),
                    values = 'black') +
  geom_boxplot(
    data = plot_data,
    aes(x = factor(x), y = value, colour = 'Hawkes'),
    outlier.shape = NA
  ) +
  scale_colour_manual(
    name = NULL,
    labels = c('Hawkes II'),
    values = c("Hawkes" = myColours[4])
  ) +
  scale_x_discrete(breaks = c(0, 25, 50 , 75)) +
  scale_y_continuous(
    breaks = c(0, 2, 4, 6),
    labels = c("0", "2", "4", "6"),
    limits = c(-2, 7)
  ) +
  labs(colour = NULL, fill = NULL) +
  xlab("t") + ylab("K(t) -  2t") +
  theme_bw() +
  theme(
    legend.position = "top",
    text = element_text(
      family = "CMU Serif",
      size = 14,
      vjust = 0.9
    )
  )

plot_data <-
  data.frame(x = seq(1, 100, 1), Hawkes = t(hawkes_arr[3, , ]))
plot_data <- melt(plot_data, id = 1)
plot_data$variable <- as.character(plot_data$variable)
plot_data$variable <- 'Hawkes'

#  Boxplot
p3 <- ggplot() +
  geom_point(aes(x = factor(1:100), y = game_values, fill = 'Observed'))  +
  scale_fill_manual(name = NULL,
                    labels = c('Observed times'),
                    values = 'black') +
  geom_boxplot(
    data = plot_data,
    aes(x = factor(x), y = value, colour = 'Hawkes'),
    outlier.shape = NA
  ) +
  scale_colour_manual(
    name = NULL,
    labels = c('Hawkes III'),
    values = c("Hawkes" = myColours[3])
  ) +
  scale_x_discrete(breaks = c(0, 25, 50 , 75)) +
  scale_y_continuous(
    breaks = c(0, 10 , 20, 30),
    labels = c("0", "10", "20", "30"),
    limits = c(-2, 30)
  ) +
  labs(colour = NULL, fill = NULL) +
  xlab("t") + ylab("K(t) -  2t") +
  theme_bw() +
  theme(
    legend.position = "top",
    text = element_text(
      family = "CMU Serif",
      size = 14,
      vjust = 0.9
    )
  )

plot_data <- data.frame(x = seq(1, 100, 1), Poisson = t(poisson_df))
plot_data <- melt(plot_data, id = 1)
plot_data$variable <- as.character(plot_data$variable)
plot_data$variable <- 'Poisson'

#  Boxplot
p4 <- ggplot() +
  geom_point(aes(x = factor(1:100), y = game_values, fill = 'Observed'))  +
  scale_fill_manual(name = NULL,
                    labels = c('Observed times'),
                    values = 'black') +
  geom_boxplot(
    data = plot_data,
    aes(x = factor(x), y = value, colour = 'Poisson'),
    outlier.shape = NA
  ) +
  scale_colour_manual(
    name = NULL,
    labels = c('Poisson'),
    values = c("Poisson" = myColours[2])
  ) +
  scale_x_discrete(breaks = c(0, 25, 50 , 75)) +
  scale_y_continuous(
    breaks = c(-2, 0, 2),
    labels = c("-  2", "0", "2"),
    limits = c(-3, 3)
  ) +
  labs(colour = NULL, fill = NULL) +
  xlab("t") + ylab("K(t) -  2t") +
  theme_bw() +
  theme(
    legend.position = "top",
    text = element_text(
      family = "CMU Serif",
      size = 14,
      vjust = 0.9
    )
  )

pdf("kfuncs.pdf", width = 10, height = 8)
p1 + p4 + p2 + p3  + plot_layout(guides = "collect") &
  theme(legend.position = 'top')
dev.off()
embed_fonts("kfuncs.pdf", outfile = "kfuncs.pdf")
