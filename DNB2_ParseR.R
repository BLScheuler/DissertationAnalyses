#DNB Step 2
#Second part of parsing


base.dir <- "../RawData/DNB/"


main_columns <- c("participant", "week", "response", "correct", "responseTime")
block_length <- 27
n_blocks <- c(8,4,4); names(n_blocks) <-  c("DualExpLoop", "AudioExpLoop", "VisualExpLoop")


all_cog <- c()
all_qualitative <- c()

for (wk in 1:3) {
  week.dir <- paste(base.dir, "Week", wk, "/", sep="")
  files <- dir(path=week.dir, pattern=".*dissDNB.*csv")
  
  for (f in files) {
    dat <- read.csv(paste(week.dir, f, sep=""))
    all_practice <- c()
    all_main <-c()
    for (loopname in c("VOFeedbackLoop", "VONoFeedbackLoop", 
                       "AOFeedbackLoop", "AONoFeedbackLoop",
                       "DualFeedbackLoop", "DualNOFeedbackLoop")) {
      signifier <- paste(loopname, ".thisN", sep="")
      practice <- subset(dat, !is.na(dat[,signifier]))[,main_columns]
      practice$trial_n <- 1:(dim(practice)[1])
      practice$block_n <- 1
      practice$trial_type <- loopname
      all_practice <- rbind(all_practice, practice)
    }
    for (loopname in c("DualExpLoop", "AudioExpLoop", "VisualExpLoop")) {
      signifier <- paste(loopname, ".thisN", sep="")
      main <- subset(dat, !is.na(dat[,signifier]))[,main_columns]
      main$trial_n <- rep(1:block_length, n_blocks[loopname])
      main$block_n <- rep(1:n_blocks[loopname], each=block_length)
      main$trial_type <- loopname
      all_main <- rbind(all_main, main)
    }
    all_cog <- rbind(all_cog, all_practice, all_main)
    
    
    qualitative_f <- subset(dat, !is.na(dat$DeviceTypeResp.started))[c("participant", "week", "age", "gender", "race.ethnicity", "cognitive.health.history", "date",
                                                                       "PhoneClick.numClicks", "TabletClick.numClicks", "MouseClick.numClicks", "LapTouchPadClick.numClicks", "TouchscreenMonitorClick.numClicks")]
    qualitative_f$device <- NA
    for (devtype in c("Phone", "Tablet", "Mouse", "LapTouchPad", "TouchscreenMonitor")) {
      if (qualitative_f[,paste(devtype, "Click.numClicks", sep="")] ==1 ) qualitative_f$device <- devtype
    }
    qualitative_f$fatigue <- subset(dat, !is.na(dat$MentalFatigue.started))[,"textbox.text"]
    all_qualitative <- rbind(all_qualitative, qualitative_f)
  }
}
all_cog$block_n <- factor(all_cog$block_n, ordered=TRUE)
all_cog$participant <- tolower(all_cog$participant)
all_cog$participant <- factor(all_cog$participant)
all_cog$response <- factor(all_cog$response)
all_cog$correct <- factor(all_cog$correct)
all_cog$trial_n <- factor(all_cog$trial_n, ordered=TRUE)
all_cog$week <- factor(all_cog$week, ordered=TRUE)
all_cog$trial_type <- factor(all_cog$trial_type)

write.csv(all_cog,file="ParseDNB.csv")