
##### Setting up the environment and data -----------------------------------------

library(readxl)
library(reshape2)
library(dplyr)
library(lubridate)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(gcplyr)

getwd() 

TestID <- read_xlsx('Data_StreamsideThermalTolerance_ThermalToleranceTrial.xlsx', sheet = 'TestID', na = 'NA') %>% filter(!is.na(StartTime))
Trial  <- read_xlsx('Data_StreamsideThermalTolerance_ThermalToleranceTrial.xlsx', sheet = 'Trial',  na = 'NA')

DAT <- merge(TestID, Trial,                         # Merges data frames 
             by = 'TestID',                         # By column Test ID
             suffixes = c('.TestID','.Trial'),      # Uses suffix to label repeated column names
             all.x = T)                             # Only keeping data for valid tests

dat <- DAT %>%
  mutate(Stream = factor(Stream, levels = c('Brushy', 'Birch', 'Littlelow', 'Poca'))) %>%
  mutate(Season = factor(Season, levels = c('Early Summer', 'Late Summer', 'Early Fall'))) %>%
  mutate(Acc_t = difftime(StartTime, TimeAtLiveWell, units = 'min'),  # Calculating acclimation time
         Rec_t = t_Rec - t_Trans) %>%                                 # Calculating recovery time
  filter(is.na(Mortality),                                            # Removes all instances of test-related mortalities
         is.na(Censor))                                               # Removes all fish that need to be censored

##### Estimating tolerance metrics for each fish in parallel ----------------------

# Set how many cores to use
n.cores <- parallel::detectCores() - 1

# Create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

# Register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

start <- Sys.time()
RES <- foreach(fish = 1:nrow(dat),
               .packages = c('dplyr', 'reshape2', 'readxl', 'lubridate', 'gcplyr'),
               .combine = 'rbind',
               .inorder = T) %dopar% {
                 if(between(dat$Date[fish], as.Date('2025-03-09'), as.Date('2025-11-02'))){
                   Heat_Tank <- read_excel(paste('ThermalToleranceTemperatureData/', dat$HT_File[fish], ' EDT.xlsx', sep = ''))[,-1]; colnames(Heat_Tank) <- c('datetime', 'Temp')
                   Reco_Tank <- read_excel(paste('ThermalToleranceTemperatureData/', dat$RT_File[fish], ' EDT.xlsx', sep = ''))[,-1]; colnames(Reco_Tank) <- c('datetime', 'Temp')
                 } else {
                   Heat_Tank <- read_excel(paste('ThermalToleranceTemperatureData/', dat$HT_File[fish], ' EST.xlsx', sep = ''))[,-1]; colnames(Heat_Tank) <- c('datetime', 'Temp')
                   Reco_Tank <- read_excel(paste('ThermalToleranceTemperatureData/', dat$RT_File[fish], ' EST.xlsx', sep = ''))[,-1]; colnames(Reco_Tank) <- c('datetime', 'Temp')
                 }
                 
                 # Creates data frame of temperatures from both heating and recovery tanks
                 Temp_dat <- merge(Heat_Tank, Reco_Tank, 
                                   by = 'datetime', 
                                   suffixes = c('Heat', 'Reco'), 
                                   all = T) 
                 
                 # Filters combined temperature data to only include time range during individual
                 # fish trials, then assigns a temperature (i.e., heating tank or recovery tank)
                 # based on the time of transfer. Finally, temperature is log-transformed and 
                 # time over the test is used as an index for next step.
                 if(!is.na(dat$t_Rec[fish])){
                   
                   Trial_dat <- Temp_dat %>%
                     filter(datetime >= dat$StartTime[fish] & datetime <= dat$t_Rec[fish]) %>%
                     rowwise() %>%
                     mutate(Temp = ifelse(datetime < dat$t_Trans[fish], yes = TempHeat, no  = TempReco),
                            log.T = log(Temp),
                            t = difftime(datetime, dat$StartTime[fish], units = 'sec'))
                 } else {
                   Trial_dat <- Temp_dat %>%
                     filter(datetime >= dat$StartTime[fish] & datetime <= dat$EndTime[fish]) %>%
                     rowwise() %>%
                     mutate(Temp = ifelse(datetime < dat$t_Trans[fish], yes = TempHeat, no  = TempReco),
                            log.T = log(Temp),
                            t = difftime(datetime, dat$StartTime[fish], units = 'sec'))
                 }
                 
                 # Find index when fish exhibited LOE (used to find CTM and censor recovery temperatures
                 # during iterative minimization procedure).
                 if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                   ind_LOE <- NA} else {
                     ind_LOE <- min(which(Trial_dat$datetime == dat$t_Trans[fish]), na.rm = T)
                   }
                 
                 # Creates blank matrix to fill with iterative estimates of area between curves 
                 # based on different starting times.
                 mat <- matrix(NA, nrow = nrow(Trial_dat), ncol = 3)
                 
                 # Fills blank matrix with estimates of area between two curves (i.e., curve 1 = 
                 # log(Temperature at time t) and the curve 2 = estimate of log(Tc))
                 for(i in 1:nrow(Trial_dat)){
                   mat[i,1] <- Trial_dat$log.T[i] - sum(Trial_dat$log.T[i:nrow(Trial_dat)]) / (nrow(Trial_dat) - i)
                 }
                 
                 mat[,2] <- abs(mat[,1])
                 
                 # Computes the moving slope of the area function across a 20 unit window of indices
                 # to verify that the area function reached a true minimum (i.e., where the sign of the derivative
                 # switches from negative to positive).
                 for(i in 11:(nrow(Trial_dat)-11)){
                   mat[i,3] <- sum((mat[(i-10):(i+10),2] - mean(mat[c((i-10):(i+10)),2])) * (c((i - 10):(i + 10)) - mean(c((i - 10):(i + 10))))) / sum((c((i - 10):(i + 10)) - mean(c((i - 10):(i + 10))))^2)
                   
                 }
                 
                 # Identifies the index of the temperature function where the area under the curve
                 # function is minimized over the range of time indices (i.e., when the temperature
                 # function reaches Tc) without influence from the endpoints or when the fish was transferred
                 # to the recovery tank.
                 if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                   ind_Tc <- NA} else {
                     extrema <- find_local_extrema(mat[1:ind_LOE,2], window_width_n = 301, return_endpoints = F, return_minima = T, return_maxima = F)
                     ind_Tc <- extrema[which(mat[extrema,2] == min(mat[extrema, 2], na.rm = T))] # Absolute minimum
                     
                   }
                 
                 # Appends the original data frame with estimates of back-transformed Tc (i.e.,
                 # when the starting temperature at which the area between the two curves equals
                 # zero), as well as estimates of CTM (time at LOE).
                 res <- matrix(NA, nrow = 1, ncol = 5)
                 
                 # Estimate of Tc
                 if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                   res[, 1] <- NA } else {
                     res[, 1] <- exp(Trial_dat$log.T[ind_Tc])
                   }
                 
                 # Estimate of CTmax
                 if(is.na(dat$t_Trans[fish])){
                   res[, 2] <- NA } else {
                     res[, 2] <- max(Trial_dat$Temp)
                   }
                 
                 # Area under the curve to make sure estimation procedure worked for each fish
                 if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                   res[, 3] <- NA } else {
                     res[, 3] <- mat[ind_Tc, 2]
                   }
                 
                 # Storing indices to evaluate problem ranges during minimization
                 if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                   res[, 4] <- NA } else {
                     res[, 4] <- ind_Tc
                   }
                 
                 # Mean temperature for the time a fish spends in the recovery tank 
                 if(is.na(dat$t_Trans[fish]) | is.na(dat$t_Rec[fish])){
                   res[, 5] <- NA } else {
                     res[, 5] <- mean(Trial_dat$Temp[Trial_dat$datetime > dat$t_Trans[fish]])
                   }
                 
                 return(res)
               }

end <- Sys.time()

# Stop Cluster
parallel::stopCluster(cl = my.cluster)

# Time Elapsed
end - start

RES <- as.data.frame(RES)
colnames(RES) <- c('Tc','CTM','AreaUnderCurve','Ind', 'MeanRecT')
dat <- cbind(dat, RES)

################################################################################
#################        Temperature logger data         #######################
################################################################################

StreamTemp_Birch     <- read_excel('StreamTemperatureLoggers/21477232 2025-10-03 11_40_36 EDT.xlsx')[,-1][-1,]; colnames(StreamTemp_Birch) <- c('datetime', 'Temp'); StreamTemp_Birch$Stream <-'Birch'
StreamTemp_Brushy    <- read_excel('StreamTemperatureLoggers/20855740 2025-09-27 14_03_27 EDT.xlsx')[,-1][-1,]; colnames(StreamTemp_Brushy) <- c('datetime', 'Temp'); StreamTemp_Brushy$Stream <- 'Brushy'
StreamTemp_Littlelow <- read_excel('StreamTemperatureLoggers/21097652 2025-10-04 13_16_44 EDT.xlsx')[,-1][-1,]; colnames(StreamTemp_Littlelow) <- c('datetime', 'Temp'); StreamTemp_Littlelow$Stream <- 'Littlelow' 
StreamTemp_Poca      <- read_excel('StreamTemperatureLoggers/20855729 2025-09-28 13_59_26 EDT.xlsx')[,-1][-1,]; colnames(StreamTemp_Poca) <- c('datetime', 'Temp'); StreamTemp_Poca$Stream <- 'Poca'

# Filling missing temperature data from 1 June 2025 to earliest collection
proxy <- read_excel('StreamTemperatureLoggers/Littlelow Proxy.xlsx')[,-1][-1,]; colnames(proxy) <- c('datetime', 'Temp')
proxy <- proxy %>%
  mutate(datetime = round_date(datetime, unit = 'hour')) %>%
  filter(as.POSIXct('2025-05-15 00:00:00') <= datetime)

pred_miss_proxy <- function(x){
  y <- merge(x %>% mutate(datetime = round_date(datetime, unit = 'hour')),
        proxy, 
        by = 'datetime',
        suffixes = c('.Measured', '.Proxy'), all = T) %>%
    mutate(Temp.Pred = predict(lm(Temp.Measured ~ Temp.Proxy, .), newdata = .)) %>%
    rowwise() %>%
    mutate(Temp = ifelse(is.na(Temp.Measured), yes = Temp.Pred, no = Temp.Measured),
           Measure = ifelse(is.na(Temp.Measured), yes = 'Predicted', no = 'Measured'))
  
  return(y)
  
}

StreamTemp_Birch     <- pred_miss_proxy(StreamTemp_Birch);     StreamTemp_Birch$Stream <- 'Birch' 
StreamTemp_Brushy    <- pred_miss_proxy(StreamTemp_Brushy);    StreamTemp_Brushy$Stream <- 'Brushy'
StreamTemp_Littlelow <- pred_miss_proxy(StreamTemp_Littlelow); StreamTemp_Littlelow$Stream <- 'Littlelow' 
StreamTemp_Poca      <- pred_miss_proxy(StreamTemp_Poca);      StreamTemp_Poca$Stream <- 'Poca'

### NOTE: Birch temp logger was likely de-watered through most of September. 

# Using Brushy as a proxy to predict the de-watered dates
pred_miss_temp <- function(Target, Proxy){
  
  X <- Target[,c('datetime','Temp', 'Stream', 'Measure')]
  Y <- Proxy[,c('datetime','Temp')]
  
  RES <- merge(X, Y, by = 'datetime', suffixes = c('.Measured', '.Proxy'), all = T) %>%
    mutate(Temp.Pred = predict(lm(Temp.Measured ~ Temp.Proxy, .), newdata = .)) %>%
    rowwise() %>%
    mutate(Temp = ifelse(is.na(Temp.Measured), yes = Temp.Pred, no = Temp.Measured),
           Measure = ifelse(Measure == 'Predicted' | is.na(Temp.Measured), yes = 'Predicted', no = 'Measured')) %>%
    dplyr::select(datetime, Temp.Measured, Stream, Temp.Proxy, Temp.Pred, Temp, Measure)
  
  return(RES)
  
}

StreamTemp_Birch$Temp[between(StreamTemp_Birch$datetime,as.POSIXct('2025-09-07'), as.POSIXct('2025-09-27'))] <- NA
StreamTemp_Birch <- pred_miss_temp(StreamTemp_Birch, StreamTemp_Brushy)
StreamTemp_Birch$Temp[StreamTemp_Birch$datetime %in% as.POSIXct(c('2025-07-02 10:00', '2025-09-06 23:00', '2025-09-07 00:00'))] <- NA

# Testing prediction accuracy with RMSE

ACCURACY <- function(x){
  X <- filter(x, Measure == 'Measured', !is.na(Temp.Pred))
  pred <- X$Temp.Pred
  actu <- X$Temp
  
  RMSE <- sqrt(sum((pred - actu) ^ 2) / length(pred))
  R2 <- summary(lm(pred ~ actu))$r.squared
  
  out <- list(RMSE = RMSE, R2 = R2)
  
  return(out)
}

ACCURACY(StreamTemp_Birch)
ACCURACY(StreamTemp_Brushy)
ACCURACY(StreamTemp_Littlelow)
ACCURACY(StreamTemp_Poca)

  
# Figure of stream temperatures ------------------------------------------------
StreamTemp_dat <- rbind(StreamTemp_Birch, StreamTemp_Brushy, StreamTemp_Littlelow, StreamTemp_Poca) 
StreamTemp_dat$Drainage  <- factor(StreamTemp_dat$Stream, levels = c('Brushy','Birch', 'Littlelow', 'Poca'), labels = c('Potomac Drainage','Ohio Drainage','Potomac Drainage','Ohio Drainage'))
StreamTemp_dat$Elevation <- factor(StreamTemp_dat$Stream, levels = c('Brushy','Birch', 'Littlelow', 'Poca'), labels = c('Low Elevation','Low Elevation','High Elevation','High Elevation'))
  
# Stream labels for each panel
annotation_data <- data.frame(
  Stream = c('Brushy','Birch','Littlelow','Poca'),
  Drainage = c('Potomac Drainage','Ohio Drainage','Potomac Drainage','Ohio Drainage'),
  Elevation = c('Low Elevation','Low Elevation','High Elevation','High Elevation'),
  Measure = 'Predicted',
  x = rep(as.POSIXct('2025-06-01'), 4),
  y = rep(21.5, 4)
)

# Points showing when sampling was done
sample_dates <- data.frame(
  Stream = rep(c('Brushy','Birch','Littlelow','Poca'),3),
  Drainage = rep(c('Potomac Drainage','Ohio Drainage','Potomac Drainage','Ohio Drainage'), 3), 
  Elevation = rep(c('Low Elevation','Low Elevation','High Elevation','High Elevation'),3),
  Measure = 'Predicted',
  x = as.POSIXct(c('2025-06-24', '2025-07-02', '2025-06-26', '2025-06-25',
                   '2025-08-08', '2025-08-11', '2025-08-29', '2025-08-05',
                   '2025-09-27', '2025-10-03', '2025-10-04', '2025-09-28')),
  y = rep(8.5, 12))

# Put whole figure together
Figure_StreamTemp <- ggplot(data = StreamTemp_dat %>% ungroup() %>% filter(datetime >= as.POSIXct('2025-06-01')),
       aes(x = datetime, y = Temp, linetype = Measure)) + 
  geom_path() + 
  ylab('Stream Temperature (\u00b0C)') + 
  geom_text(aes(x = x, y = y, label = Stream),
            data = annotation_data,
            family = 'serif', 
            size = 6, 
            hjust = 0) +
  geom_point(aes(x = x, y = y), data = sample_dates, pch = 6, size = 3) +
  geom_vline(color = 'grey35', xintercept = as.POSIXct(c('2025-07-16 00:00','2025-09-01 00:00'))) +
  coord_cartesian(ylim = c(8.5,22)) +
  facet_grid(vars(Elevation), vars(Drainage)) +
  theme_bw() +
  theme(text = element_text(family = 'serif'), 
        strip.text   = element_text(size = 14), 
        axis.text    = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text  = element_text(size = 12),
        legend.title = element_blank(),
        legend.key.width = unit(0.9, 'cm'),
        legend.position = 'bottom') 

StreamTemp_Birch <- StreamTemp_Birch[complete.cases(StreamTemp_Birch$Temp),]

### Stream temperature metrics for analysis ------------------------------------
CDH_Sum <- function(x){
  
  require(zoo)
  
  RES <- x %>%
    mutate(datetime = round_date(datetime, unit = 'hour')) %>%
    rowwise() %>%
    mutate(Over15 = ifelse(Temp > 15, yes = 1, no = 0) * (Temp - 15),
           Over16 = ifelse(Temp > 16, yes = 1, no = 0) * (Temp - 16),
           Over17 = ifelse(Temp > 17, yes = 1, no = 0) * (Temp - 17),
           Over18 = ifelse(Temp > 18, yes = 1, no = 0) * (Temp - 18),
           Over19 = ifelse(Temp > 19, yes = 1, no = 0) * (Temp - 19),
           Over20 = ifelse(Temp > 20, yes = 1, no = 0) * (Temp - 20)) %>%
    ungroup() %>%
    mutate(CDH_15 = cumsum(Over15),
           CDH_16 = cumsum(Over16),
           CDH_17 = cumsum(Over17),
           CDH_18 = cumsum(Over18),
           CDH_19 = cumsum(Over19),
           CDH_20 = cumsum(Over20)) %>%
    mutate(Mean1d  = rollmean(Temp, k = 24 *  1, fill = NA, align = 'right'),
           Mean1wk = rollmean(Temp, k = 24 *  7, fill = NA, align = 'right'),
           Mean2wk = rollmean(Temp, k = 24 * 14, fill = NA, align = 'right'),
           Mean3wk = rollmean(Temp, k = 24 * 21, fill = NA, align = 'right'),
           Mean4wk = rollmean(Temp, k = 24 * 28, fill = NA, align = 'right')) %>%
    mutate(Mov1wk_CDH_15 = rollsum(Over15, k = 24 *  7, fill = 0, align = 'right'),
           Mov1wk_CDH_16 = rollsum(Over16, k = 24 *  7, fill = 0, align = 'right'),
           Mov1wk_CDH_17 = rollsum(Over17, k = 24 *  7, fill = 0, align = 'right'),
           Mov1wk_CDH_18 = rollsum(Over18, k = 24 *  7, fill = 0, align = 'right'),
           Mov1wk_CDH_19 = rollsum(Over19, k = 24 *  7, fill = 0, align = 'right'),
           Mov1wk_CDH_20 = rollsum(Over20, k = 24 *  7, fill = 0, align = 'right')) %>% 
    mutate(Mov2wk_CDH_15 = rollsum(Over15, k = 24 * 14, fill = 0, align = 'right'),
           Mov2wk_CDH_16 = rollsum(Over16, k = 24 * 14, fill = 0, align = 'right'),
           Mov2wk_CDH_17 = rollsum(Over17, k = 24 * 14, fill = 0, align = 'right'),
           Mov2wk_CDH_18 = rollsum(Over18, k = 24 * 14, fill = 0, align = 'right'),
           Mov2wk_CDH_19 = rollsum(Over19, k = 24 * 14, fill = 0, align = 'right'),
           Mov2wk_CDH_20 = rollsum(Over20, k = 24 * 14, fill = 0, align = 'right')) %>% 
    dplyr::select(Stream, datetime, 
                  Mean1d, Mean1wk, Mean2wk, Mean3wk, Mean4wk,
                  CDH_15, CDH_16, CDH_17, CDH_18, CDH_19, CDH_20,
                  Mov1wk_CDH_15, Mov1wk_CDH_16, Mov1wk_CDH_17, Mov1wk_CDH_18, Mov1wk_CDH_19, Mov1wk_CDH_20,
                  Mov2wk_CDH_15, Mov2wk_CDH_16, Mov2wk_CDH_17, Mov2wk_CDH_18, Mov2wk_CDH_19, Mov2wk_CDH_20)

  RES
  
}

# Cumulative and moving cumulative degree hour summaries across all streams
AllStreamTemps <- mget(ls()[startsWith(ls(),'StreamTemp')][c(1,2,4,5)])
AllStreamSums  <- lapply(AllStreamTemps, FUN = CDH_Sum) %>%
  do.call('rbind', .)

# Merge with thermal tolerance data
dat <- dat %>%
  rowwise() %>%
  mutate(StartTime = if_else(TestID %in% c('CTM_Brushy_2025-09-27_3','CTM_Birch_2025-10-03_2'), # Workaround to deal with loggers not recording temperature just as tests were started
                             true = StartTime - 3600, 
                             false = StartTime)) %>% 
  mutate(datetime = round_date(StartTime, unit = 'hour')) %>%
  ungroup() %>%
  left_join(., AllStreamSums, by = c('Stream' = 'Stream', 'datetime' = 'datetime'))


################################################################################
################# Model Comparisons to assess predictors #######################
################################################################################

library(AICcmodavg)
library(tidyverse)
library(car)

extract_form <- function(x){

  modnames <- x$Modnames
  
  n = length(modnames)
  v = character(length = n)
  
  for(i in 1:n){
    
    X <- get(modnames[[i]])
    
    v[i] <- paste(X$call)[2]
  }
  
  return(v)
  
}

# Fit comparable models evaluating Tc ------------------------------------------
Tc_mod001 <- lm(Tc ~ 1              , dat)
Tc_mod002 <- lm(Tc ~ 1 + TL         , dat)
Tc_mod003 <- lm(Tc ~ 1 +      Stream, dat)
Tc_mod004 <- lm(Tc ~ 1 + TL + Stream, dat)

Tc_mod005 <- lm(Tc ~ 1               + StreamTemperature, dat)
Tc_mod006 <- lm(Tc ~ 1 + TL          + StreamTemperature, dat)
Tc_mod007 <- lm(Tc ~ 1 +      Stream + StreamTemperature, dat)
Tc_mod008 <- lm(Tc ~ 1 + TL + Stream + StreamTemperature, dat)

Tc_mod009 <- lm(Tc ~ 1               + Mean1d, dat)
Tc_mod010 <- lm(Tc ~ 1 + TL          + Mean1d, dat)
Tc_mod011 <- lm(Tc ~ 1 +      Stream + Mean1d, dat)
Tc_mod012 <- lm(Tc ~ 1 + TL + Stream + Mean1d, dat)

Tc_mod013 <- lm(Tc ~ 1               + Mean1wk, dat)
Tc_mod014 <- lm(Tc ~ 1 + TL          + Mean1wk, dat)
Tc_mod015 <- lm(Tc ~ 1 +      Stream + Mean1wk, dat)
Tc_mod016 <- lm(Tc ~ 1 + TL + Stream + Mean1wk, dat)

Tc_mod017 <- lm(Tc ~ 1               + Mean2wk, dat)
Tc_mod018 <- lm(Tc ~ 1 + TL          + Mean2wk, dat)
Tc_mod019 <- lm(Tc ~ 1 +      Stream + Mean2wk, dat)
Tc_mod020 <- lm(Tc ~ 1 + TL + Stream + Mean2wk, dat)

Tc_mod021 <- lm(Tc ~ 1               + Mean4wk, dat)
Tc_mod022 <- lm(Tc ~ 1 + TL          + Mean4wk, dat)
Tc_mod023 <- lm(Tc ~ 1 +      Stream + Mean4wk, dat)
Tc_mod024 <- lm(Tc ~ 1 + TL + Stream + Mean4wk, dat)

Tc_listmod <- mget(ls()[startsWith(ls(),'Tc_mod0')])
Tc_modcomp <- aictab(Tc_listmod, second.ord = T)
Tc_modcomp <- as.data.frame(Tc_modcomp)
Tc_modcomp <- Tc_modcomp %>% add_column(Formula = extract_form(Tc_modcomp), .after = 1)
Tc_modcomp

# Summary of top model(s)
Tc_top_mod <- get(Tc_modcomp[1,1])
summary(Tc_top_mod)
Anova(Tc_top_mod)
# TL is significant; Stream is significant; StreamTemperature is significant

# Fit comparable models evaluating CTM -----------------------------------------
CTM_mod001 <- lm(CTM ~ 1              , dat)
CTM_mod002 <- lm(CTM ~ 1 + TL         , dat)
CTM_mod003 <- lm(CTM ~ 1 +      Stream, dat)
CTM_mod004 <- lm(CTM ~ 1 + TL + Stream, dat)

CTM_mod005 <- lm(CTM ~ 1               + StreamTemperature, dat)
CTM_mod006 <- lm(CTM ~ 1 + TL          + StreamTemperature, dat)
CTM_mod007 <- lm(CTM ~ 1 +      Stream + StreamTemperature, dat)
CTM_mod008 <- lm(CTM ~ 1 + TL + Stream + StreamTemperature, dat)

CTM_mod009 <- lm(CTM ~ 1               + Mean1d, dat)
CTM_mod010 <- lm(CTM ~ 1 + TL          + Mean1d, dat)
CTM_mod011 <- lm(CTM ~ 1 +      Stream + Mean1d, dat)
CTM_mod012 <- lm(CTM ~ 1 + TL + Stream + Mean1d, dat)

CTM_mod013 <- lm(CTM ~ 1               + Mean1wk, dat)
CTM_mod014 <- lm(CTM ~ 1 + TL          + Mean1wk, dat)
CTM_mod015 <- lm(CTM ~ 1 +      Stream + Mean1wk, dat)
CTM_mod016 <- lm(CTM ~ 1 + TL + Stream + Mean1wk, dat)

CTM_mod017 <- lm(CTM ~ 1               + Mean2wk, dat)
CTM_mod018 <- lm(CTM ~ 1 + TL          + Mean2wk, dat)
CTM_mod019 <- lm(CTM ~ 1 +      Stream + Mean2wk, dat)
CTM_mod020 <- lm(CTM ~ 1 + TL + Stream + Mean2wk, dat)

CTM_mod021 <- lm(CTM ~ 1               + Mean4wk, dat)
CTM_mod022 <- lm(CTM ~ 1 + TL          + Mean4wk, dat)
CTM_mod023 <- lm(CTM ~ 1 +      Stream + Mean4wk, dat)
CTM_mod024 <- lm(CTM ~ 1 + TL + Stream + Mean4wk, dat)

CTM_listmod <- mget(ls()[startsWith(ls(),'CTM_mod0')])
CTM_modcomp <- aictab(CTM_listmod, second.ord = F)
CTM_modcomp <- as.data.frame(CTM_modcomp)
CTM_modcomp <- CTM_modcomp %>% add_column(Formula = extract_form(CTM_modcomp), .after = 1)

# Summary of top model(s)
CTM_top_mod <- get(CTM_modcomp[1,1])
summary(CTM_top_mod)
Anova(CTM_top_mod)
# TL is not significant; Mean4wk is significant

##### NOTES ON RESULTS FROM Tc & CTM ANALYSIS 

# Acclimation effects
#    Tc could be more related to short-term physiology (Time-of and mean Temp for 1 day)
#    CTM appears to be more related to long-term tolerance (numerical time-frame not as important)

# Body size relationship still exists
#    Explanations: Thermal inertia, ontogeny, heat weakening

# Compared to Stewart et al. (2024)
#    Much lower R2 values
#    Could be due to greater range of CTM values within the same assay OR higher ramping rate
#        Highlights importance of standardization among studies

##### Conducting ANCOVA --------------------------------------------------------
library(rstatix)
library(car)

# 1. Testing if interaction term is significant (homogeneity of regression slopes)
model_int_test <- lm(Tc ~ TL * StreamTemperature * Stream, dat)
Anova(model_int_test, type = 3) # type = 3 because order does not matter
# No significant interaction

# 2. Check if residuals follow normal distribution
shapiro.test(resid(Tc_top_mod))
# Significantly differs from normal distribution -> FAIL!!!
# However, due to sample size (N > 10 per group), central limit theorem applies...

# 3. Homogeneity of variances
levene_test(resid ~ as.factor(Stream), data = dat %>% mutate(resid = resid(Tc_top_mod)))
# Not significant -> PASS

# 4. Check for outliers
dat %>%
  mutate(resid = resid(Tc_top_mod)) %>%
  mutate(std.resid = scale(resid)) %>%
  filter(std.resid > 3) %>%
  as.data.frame()
# No outliers -> PASS


## ANOVA of top model to identify significant predictors
Anova(Tc_top_mod, type = 2) # Type = 2 because no interaction in model
# All predictors are significant


##### Post-hoc comparison among streams ----------------------------------------
library(emmeans)

## Goal: Identify which stream is significantly different from the rest

# Comparing streams when at means for Wt & StreamTemperature
emm.Tc_ANCOVA <- emmeans(Tc_top_mod, 
                      specs = c('Stream', 'StreamTemperature', 'TL'), 
                      at = list(StreamTemperature = mean(dat$StreamTemperature),
                                TL = mean(dat$TL)),
                      correction = 'bonferroni')

# Comparing streams when at different values for Weight when holding Stream Temperature at mean
emm.Tc_ANCOVA.TL <- emmeans(Tc_top_mod, 
                         specs = c('Stream', 'StreamTemperature', 'TL'), 
                         at = list(StreamTemperature = mean(dat$StreamTemperature),
                                   TL = seq(min(dat$TL), max(dat$TL), length.out = 50)),
                         correction = 'bonferroni')

# Comparing streams when at different values for Stream Temperature when holding Weight at mean
emm.Tc_ANCOVA.Temp <- emmeans(Tc_top_mod, specs = c('Stream', 'StreamTemperature', 'TL'), 
                           at = list(StreamTemperature = seq(min(dat$StreamTemperature), max(dat$StreamTemperature), length.out = 50),
                                     TL = mean(dat$TL)),
                           correction = 'bonferroni')

##### Plot showing significant predictors of Tc --------------------------------
library(ggpubr)

# Plot with all values possible for Stream Temperature
Tc_v_Temp <- ggplot(data = as.data.frame(emm.Tc_ANCOVA.Temp) %>%
         mutate(Stream = factor(Stream, levels = c('Brushy', 'Birch', 'Littlelow', 'Poca')),
                Tc = emmean),
       aes(x = StreamTemperature, y = Tc, fill = Stream, color = Stream)) +
  geom_point(data = dat, 
             aes(x = StreamTemperature, y = Tc, color = Stream, shape = Season),
             size = 2) + 
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = Stream), color = NA, alpha = 0.2) +
  geom_path(linewidth = 1) +
  scale_x_continuous(breaks = seq(10, 18, by = 2), minor_breaks = seq(10, 18, by = 1)) +
  scale_y_continuous(breaks = seq(20, 30, by = 2), minor_breaks = seq(20, 30, by = 1)) +
  coord_cartesian(ylim = c(20.0, 28.5)) +
  theme_bw() +
  xlab('Stream Temperature (\u00b0C)') + 
  ylab(bquote(italic(T[c])~('\u00b0C'))) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) + 
  theme(text = element_text(family = 'serif'),
        axis.text    = element_text(size = 12),
        axis.title   = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        legend.key.width = unit(1.2, "cm"))

# Plot with all values possible for Weight
Tc_v_TL <- ggplot(data = as.data.frame(emm.Tc_ANCOVA.TL) %>% 
         mutate(Stream = factor(Stream, levels = c('Brushy', 'Birch', 'Littlelow', 'Poca')),
                Tc = emmean),
       aes(x = TL, y = Tc, fill = Stream, color = Stream)) +
  geom_point(data = dat, aes(x = TL, y = Tc, color = Stream, shape = Season),
             size = 2) + 
  geom_path(linewidth = 1) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = Stream), color = NA, alpha = 0.2) +
  scale_x_continuous(breaks = seq(0, 300, by = 50), minor_breaks = seq(0, 300, by = 25)) +
  scale_y_continuous(breaks = seq(20, 30, by = 02), minor_breaks = seq(20, 30, by = 01)) +
  coord_cartesian(ylim = c(20.0, 28.5)) +
  theme_bw() +
  xlab('Total Length (mm)') + 
  ylab(bquote(italic(T[c])~('\u00b0C'))) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         shape = guide_legend(override.aes = list(size = 3))) + 
  theme(text = element_text(family = 'serif'),
        axis.text    = element_text(size = 12),
        axis.title   = element_text(size = 16),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        legend.key.width = unit(1.2, "cm"))

Figure_Tc_v_Pred <- ggarrange(Tc_v_Temp, Tc_v_TL, 
          nrow = 1, 
          common.legend = T, 
          legend = 'right',
          labels = 'AUTO',
          label.x = c(0.15, 0.07), label.y = 0.98,
          font.label = list(family = 'serif', size = 20))


##### Plot showing significant predictors of CTM -------------------------------

CTM_v_T4week <- ggplot(dat, aes(x = Mean4wk, y = CTM)) + 
  geom_point(aes(color = Stream, shape = Season), size = 2) + 
  stat_smooth(method = lm, color = 'black', fill = 'grey80', alpha = 0.5) +
  xlab(bquote(bar(T)['4week']~('\u00b0C'))) +
  ylab(bquote(CT[max]~('\u00b0C'))) +
  guides(color = guide_legend(override.aes = list(size = 3), order = 1),
         shape = guide_legend(override.aes = list(size = 3), order = 2)) + 
  theme_bw() + 
  theme(text = element_text(family = 'serif'),
        axis.text    = element_text(size = 12),
        axis.title   = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12))

CTM_v_TL <- ggplot(dat, aes(x = TL, y = CTM)) + 
  geom_point(aes(color = Stream, shape = Season), size = 2) + 
  stat_smooth(method = lm, color = 'black', fill = 'grey80', alpha = 0.5) +
  theme_bw() + 
  xlab('Total Length (mm)') +
  ylab(bquote(CT[max]~('\u00b0C'))) +
  guides(color = guide_legend(override.aes = list(size = 3), order = 1),
         shape = guide_legend(override.aes = list(size = 3), order = 2)) + 
  theme(text = element_text(family = 'serif'),
        axis.text    = element_text(size = 12),
        axis.title   = element_text(size = 16),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12))

Figure_CTM_v_Pred <- ggarrange(CTM_v_T4week, CTM_v_TL, 
                              nrow = 1, 
                              common.legend = T, 
                              legend = 'right',
                              labels = 'AUTO',
                              label.x = c(0.15, 0.07), label.y = 0.98,
                              font.label = list(family = 'serif', size = 20))


################################################################################
#####################         Ceiling analysis         #########################
################################################################################

Temp <- seq(0,35, by = 0.1)
min_TL <- min(dat$TL)
max_TL <- max(dat$TL)
TL <- seq(min_TL, max_TL, length.out = 50)

Tc_preds <- expand.grid(Int = 1,
                        TL = TL,
                        Cat1 = 0,
                        Cat2 = 0,
                        Cat3 = 0,
                        Temp = Temp)
CTM_preds <- expand.grid(Int = 1,
                         TL = TL,
                         Temp = Temp)

Tc_coef <- matrix(coef(Tc_top_mod), ncol = 1)
CTM_coef <-matrix(coef(CTM_top_mod), ncol = 1)

Tc_predicted <- cbind(Tc_preds, Tc = as.matrix(Tc_preds) %*% Tc_coef)
CTM_predicted <- cbind(CTM_preds, CTM = as.matrix(CTM_preds) %*% CTM_coef)

Tc_predicted %>% 
  as.data.frame() %>% 
  mutate(diff = Tc - Temp) %>%
  arrange(TL, Temp) %>%
  mutate(Change = ifelse(sign(diff) != sign(lag(diff)), yes = 1, no = 0)) %>%
  group_by(TL) %>%
  filter(Change == 1 & diff < 0) %>% 
  ungroup() %>%
  summarize(mean = mean(Temp), min = min(Temp), max = max(Temp))

CTM_predicted %>% 
  as.data.frame() %>% 
  mutate(diff = CTM - Temp) %>%
  arrange(TL, Temp) %>%
  mutate(Change = ifelse(sign(diff) != sign(lag(diff)), yes = 1, no = 0)) %>%
  group_by(TL) %>%
  filter(Change == 1 & diff < 0) %>% 
  ungroup() %>%
  summarize(mean = mean(Temp), min = min(Temp), max = max(Temp))

new_labels <- as_labeller(c(Tc = "italic(T[c])", CTM = "italic(CT[max])"), default = label_parsed)

Ceiling_Plot <- ggplot(data = rbind(filter(Tc_predicted, TL %in% c(min_TL, TL[26], max_TL)) %>%
                      dcast(Temp ~ TL, value.var = 'Tc') %>%
                      `colnames<-`(c('Temp', "LCL", "Mean", "UCL")) %>%
                      mutate(Variable = 'Tc'),
                    filter(CTM_predicted, TL %in% c(min_TL, TL[26], max_TL)) %>%
                      dcast(Temp ~ TL, value.var = 'CTM') %>%
                      `colnames<-`(c('Temp', "LCL", "Mean", "UCL")) %>%
                      mutate(Variable = 'CTM')) %>% 
                      mutate(Variable = factor(Variable, levels = c("Tc", 'CTM'))),
       aes(x = Temp)) + 
  coord_cartesian(xlim = c(0,35), ylim = c(0,35)) +
  geom_path(aes(y = Mean)) + 
  geom_path(aes(y = LCL), linetype = 2) + 
  geom_path(aes(y = UCL), linetype = 2) + 
  geom_abline(intercept = 0, slope = 1) +  
  facet_wrap(~Variable, labeller = as_labeller(new_labels),) +
  xlab('Acclimation Temperature Metric (\u00b0C)') + 
  ylab('Thermal Tolerance Metric (\u00b0C)') +
  theme_bw() + 
  theme(text = element_text(family = 'serif'),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16))

ggsave('Results/CeilingPlot.jpg', plot = Ceiling_Plot, width = 8, height = 6, units = 'in', dpi = 'print')
       
################################################################################
######################         Important notes         #########################
################################################################################

### Results

# 1) Stream temperature and body weight were significant predictors of thermal 
#    tolerance in field trials.
# 2) Brushy (lowest elevation stream) has narrowly significantly different thermal 
#    tolerance than other streams when controlling for stream temperature and 
#    individual weights.

### Interpretation and impact

# 1) Adds to growing body of literature demonstrating that acclimation temperature 
#    is an important predictor for thermal tolerance of ectotherms.
# 2) Suggests that older fish that have likely experienced thermal stressors in 
#    past may be less tolerant in the future (did not assign individual ages). 
#    There were no body size differences in previous study when ages were all the 
#    same but DID decrease with increasing number of thermal challenges... (DO NOT
#    MAKE TOO STRONG OF A CLAIM HERE!!!)
# 3) Brushy's difference is likely a function of small range of stream temperatures
#    in analyses (next top model didn't include Stream), and is not interesting. 
# 4) Acclimation effects on Tc will eventually fail to keep up with negative 
#    effects of stream temperature (i.e., Stream temperature will exceed Tc curve)


################################################################################
#######################         Export Results         #########################
################################################################################
library(writexl)

# Exporting model selection and top model summaries
write_xlsx(x = list('TcModelComparisons' = Tc_modcomp,
                    'TcANOVAResults' = as.data.frame(Anova(Tc_top_mod) %>% mutate(Parameter = rownames(.), .before = 1)),
                    'TcTopModelSummary'  = tidy(Tc_top_mod),
                    'CTMModelComparisons' = CTM_modcomp,
                    'CTMANOVAResults' = as.data.frame(Anova(CTM_top_mod) %>% mutate(Parameter = rownames(.), .before = 1)),
                    'CTMTopModelSummary'  = tidy(CTM_top_mod)),
           path = 'Results/ModelSummaries_Export.xlsx')

# Exporting figures
ggsave('Results/StreamTemperature.jpeg', 
       Figure_StreamTemp, 
       width = 8, 
       height = 7, 
       unit = 'in',
       dpi= 'print')

ggsave('Results/Tc v. Predictors.jpeg', 
       Figure_Tc_v_Pred, 
       width = 11, 
       height = 8, 
       unit = 'in',
       dpi= 'print')

ggsave('Results/CTM v. Predictors.jpeg', 
       Figure_CTM_v_Pred, 
       width = 11, 
       height = 8, 
       unit = 'in',
       dpi= 'print')

