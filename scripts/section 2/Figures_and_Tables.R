## This script generates Tables 1,2,3 and Figures 1,2
## in sections 1 and 2 of the paper

# Packages
library(data.table)
library(xtable)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)
library(ggpubr)
library(extrafont)

# Load the data
load_data <- function(path) { 
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)
  
  match.ids <- vapply(files,function(x){
    as.numeric(unlist(strsplit(x, "[_.]"))[2])
  },FUN.VALUE = 0)
  tables <- lapply(files, read.csv)
  
  for(i in 1:length(tables)){
    x <- data.frame(tables[i])
    x$match_id <- match.ids[i]
    tables[[i]] <- x
  }

  dplyr::bind_rows(tables)
}

dir_path = '/Users/santhoshnarayanan/Documents/PhD/EngPr/'

# Meta Data
team.df <- as.data.table(
  read.csv(paste0(dir_path,"fixture_info_engpr.csv"))[,c(3,7,8,10,11)])

# Data Import
data <-
  load_data(paste0(dir_path,"season1"))

# Basic cleaning

# Remove events which do not have associated location
data <-  data[-which(data$type %in% c("PenaltyFaced","ChanceMissed","BlockedPass", "ShieldBallOpp", 
                                      "FormationChange", "FormationSet", "KeeperSweeper", "Card", 
                                      "Start", "SubstitutionOff", "SubstitutionOn", "End","OffsideGiven","Error")),]

# Remove events which come in pairs (duplicates)
data <-  data[-which(data$type %in% c("CornerAwarded","Foul","Aerial","Challenge") & data$outcome == "Unsuccessful"),]
data <-  data[-which(data$type == "OffsideProvoked"),]
data$second <- as.integer(data$second)
data$player_id <- as.factor(data$player_id)

# Table 1 - data set snapshot
team.df[team.df$fixture_id == 1483492,]
# Southampton vs West Ham United on 15/9/2013
snapData <- head(data[data$match_id == 1483492, 
                      c(9, 5, 10, 8, 11, 12:13, 6, 2:3)], n=9)
snapData[snapData$type != 'Pass', 'outcome'] <- ''
snapData[snapData$type != 'Pass', 'end_x'] <- ''
snapData[snapData$type != 'Pass', 'end_y'] <- ''
snapData$team_id[snapData$team_id == 670] <- as.integer(14)
snapData$team_id[snapData$team_id == 684] <- as.integer(20)
print(xtable(snapData, digits = 1), include.rownames = F)

# Table 2 - Season 1 Teams
teams_s1 <- data.table(data.frame( team_id = 1:20, team_name = sort(unique(team.df$team1[1:380])) ) )
print(xtable(cbind(teams_s1[1:10,],teams_s1[11:20,])), include.rownames = F)

dataset <- data

# Table 3 - Event frequency
agg.df = aggregate.data.frame(x = data.frame(Frequency = dataset$type),by= list(Event = dataset$type),FUN = length)
agg.df = agg.df[order(agg.df$Frequency, decreasing = T),]
agg.df <- agg.df[-which(agg.df$Event == "GoodSkill"),]
row.names(agg.df) <- 1:nrow(agg.df) 
print(xtable(cbind(agg.df[1:11,],agg.df[12:22,]), digits = 1), include.rownames = F)

# data frame with pitch markers
pitch <- data.frame(x1 = c(0,0,0,100,50,0,16.5,16.5,100,83.5,83.5),
                    x2 = c(100,0,100,100,50,16.5,16.5,0,83.5,83.5,100),
                    y1 = c(0,0,100,100,0,22,22,78,22,22,78),
                    y2 = c(0,100,100,0,100,22,78,78,22,78,78))


# VISUALIZATIONS

# Select the game with events leading to Goal of the season
dummy = dataset[dataset$match_id=='1483527',]
dummy = dummy[316:328,]
dummy$time <-  60*dummy$expanded_minute + dummy$second
df = data.frame(x= NA, y = NA, team = NA, type=NA, time=NA)

# Create a new data frame with a new row between every recorded events
for(i in 1:nrow(dummy)){
  
  if(i > 1){
    if( dummy[i,"team_id"] != dummy[i-1,"team_id"]){
      row = data.frame(NA,NA,dummy[i,"team_id"], NA, NA)
      colnames(row) <- colnames(df)
      df = rbind(df,row)
    }
  }
  
  row1 = dummy[i,c("x","y","team_id","type","time")]
  row2 = dummy[i,c("end_x","end_y","team_id","type","time")]
  row2$type <- 'Touch'
  
  colnames(row1) <- colnames(row2) <- colnames(df)
  df = rbind(df,row1,row2)
  
}

# Some tweaks before plot
df = df[-1,]
df <- df[!duplicated(df[,1:3]),]
df$type[df$type != 'Pass'] <- 'Touch'
df = rbind(df,c(100, 53, 660, "End", 1026))
df = df[-c(2,11,13,15,17,19),]
df$x <- as.numeric(df$x)
df$y <- as.numeric(df$y)
df$team = factor(df$team)
df$time <- as.numeric(df$time)
df$time <- df$time - min(df$time)

amin <- 0.2
highcol <- 'red'
lowcol.hex <- as.hexmode(round(col2rgb(highcol) * amin + 255 * (1 - amin)))
lowcol <- paste0("#",   sep = "",
                 paste(format(lowcol.hex, width = 2), collapse = ""))

# Figure 1 - Goal of the season
p <- ggplot(df, aes(x = x, y = y)) + 
  geom_path(colour="blue") + 
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = pitch,colour="grey") +  
  geom_point(aes(alpha = time), colour = highcol, size = 4) +
  geom_point(aes(colour = time), alpha = 0) + 
  scale_colour_gradient(high = highcol, low = lowcol) + 
  guides(alpha = F) + 
  labs(colour = "time (s)") +
  theme_bw() + 
  scale_y_continuous(expand = c(0,0), breaks=NULL) + 
  scale_x_continuous(expand = c(0,0), breaks=NULL) + 
  theme(legend.position="bottom", axis.title= element_blank(),
        text=element_text(family="CMU Serif", size=14, vjust = 0.9))

ggsave("figures/goalpath.pdf", plot=p,  width=8, height=5)
embed_fonts("figures/goalpath.pdf", outfile="figures/goalpath.pdf")

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
 
# assignInNamespace("colour_ramp", function(colors, na.color = NA, alpha = TRUE){
#   if (length(colors) == 0) {
#     stop("Must provide at least one color to create a color ramp")
#   }
#   colorMatrix <- grDevices::col2rgb(colors, alpha = alpha)
#   structure(function(x) {
#     scales:::doColorRamp(colorMatrix, x, alpha, ifelse(is.na(na.color), 
#                                                        "", na.color))
#   }, safe_palette_func = TRUE)
# }, "scales")

my.labs <-  function(breaks){
  labels = c('low', rep('', length(breaks) - 2) ,'high')
}

my.breaks <-  function(limits){
  seq(limits[1], limits[2], length.out = 10)
}

my.heatmap <- function(heat.data, my.title = ''){
  df <- data.frame(x=heat.data$x, y=heat.data$y)
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  my.colours <- alpha(jet.colors(10), seq(0.5, 1, length.out = 10))
  
  ggplot(df, aes(x, y)) + 
    stat_density2d(aes(fill=..density..), geom="raster", contour=FALSE, n=200) +  
    scale_fill_gradientn(name = 'density  ', labels = my.labs, breaks = my.breaks, colours = my.colours, trans = 'identity') + 
    guides(alpha = F, fill = guide_colourbar(ticks = F)) +
    theme_bw() + scale_y_continuous(breaks=NULL, expand = c(0,0)) + scale_x_continuous(breaks=NULL, expand = c(0,0)) + 
    theme(axis.title= element_blank(), legend.position = 'bottom') +
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = pitch, colour="black", alpha = 0.75) + 
    ggtitle(my.title) +
    theme(plot.title = element_text(hjust = 0.5))
  
}

teams = c(660, 661, 662, 663, 670)
# Arsenal, Chelsea, Manchester United, Liverpool

p1 <- my.heatmap(data[which(data$team_id %in% teams[1] & data$match_id %in% team.df$fixture_id[team.df$team1_id == teams[1]]),])
p2 <- my.heatmap(data[which(data$team_id %in% teams[1] & data$match_id %in% team.df$fixture_id[team.df$team2_id == teams[1]]),])
p3 <- my.heatmap(data[which(data$team_id %in% teams[2] & data$match_id %in% team.df$fixture_id[team.df$team1_id == teams[2]]),])
p4 <- my.heatmap(data[which(data$team_id %in% teams[2] & data$match_id %in% team.df$fixture_id[team.df$team2_id == teams[2]]),])

row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Arsenal", angle = 90, family="CMU Serif", size=6) + theme_void()
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Chelsea", angle = 90, family="CMU Serif", size=6) + theme_void()
col1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Home", family="CMU Serif", size=6) + theme_void() 
col2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label="Away", family="CMU Serif", size=6) + theme_void() 

plotlist <- list(first = plot_spacer(), a = row1, b = row2, c = col1, d = col2, e= p1, f=p2, g=p3, h=p4)

layoutplot <- c(
  area(1, 1),
  area(2, 1, 2),
  area(3, 1, 3),
  area(1, 2),
  area(1, 3),
  area(2, 2),
  area(2, 3),
  area(3, 2),
  area(3, 3)
)

# Figure 2 - Heat maps
p <- wrap_plots(plotlist, guides = 'collect', design = layoutplot, widths = c(0.4, 4, 4), heights  = c(0.4, 4, 4)) & 
  theme(legend.position = 'bottom') & 
  theme(text=element_text(family="CMU Serif", size=14, vjust = 0.8)) 

ggsave("figures/homevsaway.pdf", plot=p,  width=10, height=8)
embed_fonts("figures/homevsaway.pdf", outfile="figures/homevsaway.pdf")

# Figure 3 - Zone
p <- ggplot() + 
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = pitch, 
               colour="grey", alpha = 0.9, linetype = 'dotted') +  
  geom_vline(xintercept = c(100/3, 200/3)) + 
  geom_text( aes(x = c(50/3, 50, 100 - 50/3), y = c(50, 50, 50), label = 1:3), size = 10) + 
  theme_bw() + 
  scale_y_continuous(expand = c(0,0), breaks=NULL) + 
  scale_x_continuous(expand = c(0,0), breaks=NULL) + 
  xlab(expression('Home team attacks' %->% '')) + ylab("") + ggtitle('Zones') + 
  theme(legend.position="bottom", plot.title = element_text(hjust = 0.5), 
        text=element_text(family="CMU Serif", size=14, vjust = 0.9)) 

ggsave("figures/zones.pdf", plot=p,  width=5, height=3)
embed_fonts("figures/zones.pdf", outfile="figures/zones.pdf")
