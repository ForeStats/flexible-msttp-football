# Data preparation script

# Load packages
library(dplyr)
library(data.table)
library(gtools)
library(xtable)

# Data load function
load_data <- function(path) {
  files <- dir(path, pattern = '\\.csv', full.names = TRUE)

  match.ids <- vapply(files, function(x) {
    as.numeric(unlist(strsplit(x, "[_.]"))[2], stringsAsFactors = TRUE)
  }, FUN.VALUE = 0)

  tables <- lapply(files, FUN = function(z) {
    df <- as.data.table(read.csv(z, stringsAsFactors = TRUE), stringsAsFactors = TRUE)
    changeCols <- c("period", "type")
    df[, (changeCols) := lapply(.SD, as.character),.SDcols = changeCols]
  })
  names(tables) <- match.ids
  dplyr::bind_rows(tables, .id = "match_id")
}

dir_path = '~/Documents/PhD/EngPr/'

# Meta Data
df <- as.data.table(
    read.csv(paste0(dir_path,"fixture_info_engpr.csv"), stringsAsFactors = TRUE)[,c(3,7,8,10,11)],
    stringsAsFactors = TRUE)

colnames(df) <- c("match_id","HT","Home","AT","Away")
df <- data.table(rbind(data.frame(match_id = df$match_id, name = df$HT,
                                  team = 'Home', team_id = df$Home, opp.name = df$AT, stringsAsFactors = TRUE),
                       data.frame(match_id = df$match_id, name = df$AT,
                                  team = 'Away', team_id = df$Away, opp.name = df$HT, stringsAsFactors = TRUE)),
                 stringsAsFactors = TRUE)
setkey(df, match_id)
changeCols <- c("match_id","team_id")
df[, (changeCols) := lapply(.SD, as.factor),.SDcols = changeCols]

# Data Import
data1 <-
  load_data(paste0(dir_path,"season1"))
data2 <-
  load_data(paste0(dir_path,"season2"))

# Formatting to factors
data <- rbind(data1,data2)
changeCols <- c("outcome","period","player_id","team_id","type","match_id")
data[, (changeCols) := lapply(.SD, as.factor),.SDcols = changeCols]
# create time field
data[, time := 60*expanded_minute + second]
# Remove offside given events without proper time field
data <- data[is.finite(time)]
# Remove events without proper location field
data <- data[!period %in% c('PreMatch','PostGame')]
data <- data[!type %in% c("SubstitutionOff",'SubstitutionOn','BlockedPass','Error','KeeperSweeper',
                          'Card','OffsideProvoked', 'FormationSet','FormationChange','ChanceMissed',
                          'GoodSkill','ShieldBallOpp','PenaltyFaced','BallRecovery', 'Aerial')]
# Remove events that come in pairs
data <- data[-which(type %in% c('Foul','CornerAwarded','Challenge') & outcome == 'Unsuccessful')]
# Merge with meta data
data <- droplevels(data)
data[, id := seq_len(.N)]

data <- data[df, on =c("match_id","team_id")]
cols <- c("id","X","match_id","period","name","opp.name","time",'team',"type","outcome","x","y","end_x","end_y")
data <- data[,.SD,.SDcols = cols]
setkey(data,id)

# Own goals
ogs <- which(data$x < 20 & data$type == 'Goal')
data[ogs, ':='(team = ifelse(team == 'Home', 'Away', 'Home'), x = 100 - x, y = 100 - y,
               end_x = 100 - end_x, end_y = 100 - end_y, name = opp.name, opp.name = name)]

for(i in 1:4){
  data[, delta := c(0,diff(time)), by = c('match_id','period')]
  data[, time := cumsum(delta), by = c('match_id','period')]

  data[, ptype := shift(type,1,type='lag'), by = c('match_id','period')]
  data[, ntype := shift(type,1,type='lead'), by = c('match_id','period')]
  data[, pdelta := shift(delta,1,type='lag'), by = c('match_id','period')]
  data[, ndelta := shift(delta,1,type='lead'), by = c('match_id','period')]
  data[, pperiod := shift(period,1,type='lag'), by = c('match_id','period')]
  data[, pteam := shift(team,1,type='lag'), by = c('match_id','period')]
  data[, nteam := shift(team,1,type='lead'), by = c('match_id','period')]
  data[, pout := shift(outcome,1,type='lag'), by = c('match_id','period')]

  # wrong terminal events
  wrong_goals <- which(data$ptype == 'Goal' & data$delta < 7)
  wrong_throws <- which(data$y %in% c(0,100) & !data$x %in% c(0,100) & data$type == "Pass"
                        & (data$delta == 0 | (data$delta < 3 & data$pdelta > 7)) )
  wrong_fks <- which(data$ptype == 'Foul' & data$delta < 3)
  wrong_corners <- which(data$ptype == 'CornerAwarded' & data$delta < 3)
  wrong_offsides <- which(data$ptype == 'OffsidePass' & data$delta < 3)
  wrong_ends <- which(data$ptype == 'End' & !data$type %in% c('Start','End')
                      & data$delta >= 0 & data$period == data$pperiod)

  # wrong same time events
  wrong_saves <- which(data$type == 'SavedShot' & data$ptype == 'Save' & data$delta == 0)
  wrong_tackles <- which(data$type %in% c('TakeOn','Dispossessed') & data$team != data$pteam
                         & data$ptype == 'Tackle' & data$delta == 0)
  wrong_clear <- which(data$ptype == 'Clearance'
                       & data$type == 'Pass' & data$outcome == 'Unsuccessful'
                       & data$team != data$pteam & data$delta == 0)
  wrong_inter <- which(data$ptype == 'Interception' & data$pout == 'Successful'
                       & data$type == 'Pass' & data$outcome == 'Unsuccessful'
                       & data$team != data$pteam & data$delta == 0)
  wrong_touch <- which(data$ptype == 'Pass' & data$pout == 'Successful'
                       & data$type == 'BallTouch'
                       & data$team != data$pteam & data$delta == 0)

  mistakes <- unique(c(wrong_throws, wrong_goals, wrong_fks, wrong_corners, wrong_ends, wrong_offsides,
                       wrong_saves, wrong_tackles, wrong_clear, wrong_inter, wrong_touch))

  data[mistakes,  ':='(id = as.integer(id - 1), time = time - delta)]
  data[mistakes - 1, id := as.integer(id + 1)]

  wrong_starts <- which(data$type == 'Start' & data$ptype == 'Pass' & data$ntype == 'Start')
  data[wrong_starts - 1, id := as.integer(id + 2)]
  data[wrong_starts + 1, id := as.integer(id - 2)]

  setkey(data,id)
  data[, delta := c(0,diff(time)), by = c('match_id','period')]
  data[, time := cumsum(delta), by = c('match_id','period')]
}

# Add Terminal Marks
data[, ptype := shift(type,1,type='lag'), by = c('match_id','period')]
# Throw-ins
throws <- which(data$y %in% c(0,100) & !data$x %in% c(0,100) & data$type == "Pass")
throwData <- data[throws]
throwData[, ':='(id = id - 0.5, end_x = x, end_y = y, type = 'Out_Throw',
                 outcome = 'Successful', time = time - delta, delta = 0)]
# goal kicks
gks <- which(data$x >= 3 & data$x <= 6 & data$y >= 36 & data$y <= 63 & data$delta > 7 & data$type == "Pass" & !data$ptype %in% c('Foul','OffsidePass'))
gkData <- data[gks]
gkData[, ':='(id = id - 0.5, end_x = x, end_y = y, type = 'Out_GK',
              outcome = 'Successful', time = time - delta, delta = 0)]

# Keeper pickups
data[, ntype := shift(type,1,type='lead'), by = c('match_id','period')]
kpks <- which( ((data$end_x >= 84 & data$outcome == 'Unsuccessful') | (data$end_x <= 16 & data$outcome == 'Successful')) &
                 data$ndelta > 9 & data$type == "Pass" & data$ntype == "Pass" & !data$end_y %in% c(0, 100) & !data$end_x %in% c(0, 100) )
kpData <- data[kpks]
kpData[outcome == 'Unsuccessful', ':='(team = ifelse(team == 'Home', 'Away', 'Home'), x = 100 - x, y = 100 - y,
                                       end_x = 100 - end_x, end_y = 100 - end_y, name = opp.name, opp.name = name)]
kpData[, ':='(id = id + 0.5, x = end_x, y = end_y, type = 'KeeperPickup',
              outcome = 'Successful', time = time + 1, delta = 1)]

data = rbind(data, throwData, gkData, kpData)
data[, ptype := NULL]
setkey(data,id)
data[, id := seq_len(.N)]
setkey(data,id)

# Away location swap
data[team == 'Away', ':='(x = 100 - x, y = 100 - y,
                          end_x = 100 - end_x, end_y = 100 - end_y)]

# Terminal checks
data[, delta := c(0,diff(time)), by = c('match_id','period')]
data[, time := cumsum(delta), by = c('match_id','period')]
data[, ndelta := shift(delta,1,type='lead'), by = c('match_id','period')]
data[, ntype := shift(type,1,type='lead'), by = c('match_id','period')]

wrong_throws = which(data$type == 'Out_Throw' & data$ndelta < 3)
wrong_corners = which(data$type == 'CornerAwarded' & (data$ntype != 'Pass' | data$ndelta < 4) )
wrong_fouls = which(data$type == 'Foul' &
                      (!data$ntype %in% c('Pass','Goal','SavedShot','MissedShots','ShotOnPost','OffsidePass')
                       | data$ndelta < 4) )
wrong_offside = which(data$type == 'OffsidePass' & (data$ntype != 'Pass' | data$ndelta < 4) )
wrong_gks = which(data$type == 'Out_GK' & (data$ntype != 'Pass' | data$ndelta < 4) )
wrong_goals = which(data$type == 'Goal' & (data$ntype != 'Pass' | data$ndelta < 4) )

wrongs_others = which(!data$type %in% c('Out_Throw', 'CornerAwarded', 'Foul', 'Claim',
                                        'OffsidePass', 'Out_GK', 'Goal', 'KeeperPickup') &
                        data$ndelta > 9)

data[wrongs_others + 1, delta := round(runif(length(wrongs_others), 5, 10), 0)]

long_gaps <- which(data$delta > 60)
data[delta > 60, delta := round(runif(length(long_gaps), 30, 60), 0)]
data[, time := cumsum(delta), by = c('match_id','period')]

data[, ':='(pdelta = NULL, ndelta = NULL, ntype = NULL,
            pperiod = NULL, pteam = NULL, nteam = NULL, pout = NULL)]

# Set Marks
data[type == 'Out_Throw', mark := 'Out_Throw']
data[type == 'Out_GK', mark := 'Out_GK']
data[type == 'CornerAwarded', mark := 'Out_Corner']
data[type == 'Foul', mark := 'Foul']
data[type == 'Goal', mark := 'Goal']
data[type == 'OffsidePass', ':='(mark = 'Pass', outcome = 'Offside')]
data[type == 'End', mark := 'End']

# Removing starts and double ends
ends = which(data$type == 'End')
data = data[-ends[even(ends)]]
starts = which(data$type == 'Start')
data = data[-starts]

# Handling simultaneous event times
data[, delta := c(0,diff(time)), by = .(match_id, period)]
data[delta == 0, delta := 1]
data[, time := cumsum(delta), by = .(match_id, period)]
data[, delta := NULL]

# Assign marks
data[type %in% c('Clearance','Punch'), mark := 'Clear']
data[type %in% c('SavedShot','ShotOnPost','MissedShots'), mark := 'Shot']
data[type == 'Pass', mark := 'Pass']
data[type %in% c('Interception','Smother'), mark := 'Win']
data[type == 'Save', mark := 'Save']
data[type %in% c('BallTouch','Dispossessed','CrossNotClaimed'), mark := 'Lose']
data[type %in% c('Tackle','TakeOn','Claim') & outcome == 'Unsuccessful', mark := 'Lose']
data[type == 'KeeperPickup', mark := 'Keeper']
data[type == 'Claim' & outcome == 'Successful', mark := 'Keeper']
data[type == 'TakeOn' & outcome == 'Successful', mark := 'Dribble']
data[type == 'Tackle' & outcome == 'Successful', mark := 'Win']

data[, delta := c(0,diff(time)), by = .(match_id, period)]
data[, p_mark := shift(mark, 1, type = 'lag'), by = c('match_id', 'period')]
data[, p_outcome := shift(outcome, 1, type = 'lag'), by = c('match_id', 'period')]

data[ (p_mark %in% c('Foul', 'Out_Throw', 'Out_GK', 'Out_Corner', 'Keeper', 'Goal')) |
        (p_mark == 'Pass' & p_outcome == 'Offside'), delta := 1]

data[delta > 10, delta := 10]
data[, time := cumsum(delta), by = .(match_id, period)]
table(data$delta)
data[, delta := NULL]
data[, p_mark := NULL]
data[, p_outcome := NULL]

# final data
df[, season := ifelse(substr(match_id,2,2) == 4, 1, 2)]
df = df[team == 'Home']
cols = c('season','match_id', 'name','opp.name')
df <- df[, .SD, .SDcols = cols ]
colnames(df) <- c('season','match_id','home','away')

# Some final mods
data <- data[mark != 'End']
data[, id := seq_len(.N)]

cols <- c('id',"match_id", 'period', "team","time",'mark',"outcome","x","y","end_x","end_y")
data <- data[,.SD,.SDcols = cols]
setkey(data,id)

# Post processing
data[mark == 'Pass', composite_mark := paste(team, mark, substr(outcome, 1, 1), sep='_')]
data[mark != 'Pass', composite_mark := paste(team, mark, sep='_')]

mark_levels <- c('Home_Win', 'Home_Dribble', 'Home_Pass_S', 'Home_Pass_U', 'Home_Shot', 'Home_Keeper', 'Home_Save', 'Home_Clear', 'Home_Lose',
                 'Home_Goal', 'Home_Foul', 'Home_Out_Throw', 'Home_Out_GK', 'Home_Out_Corner', 'Home_Pass_O',
                 'Away_Win', 'Away_Dribble', 'Away_Pass_S', 'Away_Pass_U', 'Away_Shot', 'Away_Keeper', 'Away_Save', 'Away_Clear', 'Away_Lose',
                 'Away_Goal', 'Away_Foul', 'Away_Out_Throw', 'Away_Out_GK', 'Away_Out_Corner', 'Away_Pass_O')
mark_labels <- 1:30

data[, mark_id := factor(composite_mark, levels = mark_levels, labels = mark_labels)]
data[, period := factor(period, labels = 1:2)]
data[, team := factor(team, labels = 1:2)]

### Start New stuff

#### Zones
data[, zone_s := ifelse(x > 66.6, 3, ifelse(x < 33.3 , 1, 2))]
data[, zone_f := ifelse(end_x > 66.6, 3, ifelse(end_x < 33.3 , 1, 2))]

#### Score
data[, score := 0]
data[, score := ifelse(mark_id == 10, 1, ifelse(mark_id == 25, -1, 0))]
data[, score := cumsum(score), by = 'match_id']
data[, score := ifelse(score < -1, -2, ifelse(score > 1, 2, score))]

#### History
data[, history := 1]
data[, history_flag := 0]
data[ (mark %in% c('Foul', 'Out_Throw', 'Out_GK', 'Out_Corner', 'Keeper', 'Goal')) |
        (mark == 'Pass' & outcome == 'Offside'), history_flag := 1]
data[, history_flag := 1+ cumsum(history_flag), by =  c('match_id', 'period')]
data[, history := cumsum(history), by =  c('match_id', 'period', 'history_flag')]
data[, history := shift(history, 1, type = 'lag', fill = 0), by = c('match_id', 'period')]
data[history > 5, history := 5]
data[, history_flag := NULL]

# Table 4
my.table = data[, .N, by = c('mark_id','composite_mark')]
setkeyv(my.table, c('mark_id','composite_mark','N'))
print(xtable(my.table), include.rownames = F)

df_season_1 <- df[season == 1]
df_season_1 <- droplevels(df_season_1)
df_season_1[, game := factor(match_id, labels = 1:380)]
df_season_1[, h_id := factor(home, labels = 1:20)]
df_season_1[, a_id := factor(away, labels = 1:20)]

df_season_2 <- df[season == 2]
df_season_2 <- droplevels(df_season_2)
df_season_2[, game := factor(match_id, labels = 1:380)]
df_season_2[, h_id := factor(home, labels = 1:20)]
df_season_2[, a_id := factor(away, labels = 1:20)]

data_season_1 <- data[df_season_1, on = 'match_id']
data_season_2 <- data[df_season_2, on = 'match_id']

cols <- c("game", 'period', 'time', 'mark_id', 'zone_s', 'zone_f', 'score', 'history')
data_season_1 <- data_season_1[,.SD,.SDcols = cols]
setkeyv(data_season_1, c('game', 'period', 'time'))

data_season_2 <- data_season_2[,.SD,.SDcols = cols]
setkeyv(data_season_2, c('game', 'period', 'time'))

colnames(data_season_1)[4] <- colnames(data_season_2)[4] <- 'mark'
data_season_1$mark = as.numeric(data_season_1$mark)
data_season_2$mark = as.numeric(data_season_2$mark)

list = c('data', 'df', 'cols')

save(list = c('data_season_1', 'data_season_2', 'df_season_1', 'df_season_2', 'mark_levels'),
     file = 'matchData.RData')
