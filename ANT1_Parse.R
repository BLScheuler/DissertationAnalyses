#ANT Step 1
#Parsing

#Set working directory
#Pull in raw data

prac_columns <- c("participant", "week", "cue", "tar", "corrAns", "targOrientation", "pracLoop.thisN", "resp.rt", "resp.corr")
main_columns <- c("participant", "week", "cue", "tar", "corrAns", "targOrientation", "blockLoop.thisRepN", "mainLoop.thisN", "resp.rt", "resp.corr")

all_cog <- c()
all_qualitative <- c()

for (wk in 1:3) {
  week.dir <- paste(base.dir, "Week", wk, "/", sep="")
  files <- dir(path=week.dir, pattern=".*dissANT.*csv")
  
  for (f in files) {
    dat <- read.csv(paste(week.dir, f, sep=""))
    practice <- subset(dat, !is.na(pracLoop.thisTrialN))[,prac_columns]
    practice$trial_type <- "practice"
    practice$block_n <- -1
    names(practice)[grep("thisN", names(practice))] <- "trial_n"
    all_cog <- rbind(all_cog, practice)
    main <- subset(dat, !is.na(mainLoop.thisTrialN))[,main_columns]
    main$trial_type <- "main"
    names(main)[grep("thisN", names(main))] <- "trial_n"
    names(main)[grep("thisRepN", names(main))] <- "block_n"
    all_cog <- rbind(all_cog, main)
    
    qualitative_f <- subset(dat, !is.na(dat$FatigueQ.start))[c("participant", "week", "age", "gender", "race.ethnicity", "cognitive.health.history", "date","textbox.text")]
    all_qualitative <- rbind(all_qualitative, qualitative_f)
    
  }
}
all_cog$block_n <- factor(all_cog$block_n + 1, ordered=TRUE)
all_cog$participant <- tolower(all_cog$participant)
all_cog$participant <- factor(all_cog$participant)
all_cog$week <- factor(all_cog$week, ordered=TRUE)
all_cog$cue <- factor(all_cog$cue)
all_cog$tar <- factor(all_cog$tar)
all_cog$corrAns <- factor(all_cog$corrAns)
all_cog$trial_n <- factor(all_cog$trial_n, ordered=TRUE)
all_cog$targOrientation <- factor(all_cog$targOrientation)
all_cog$trial_type <- factor(all_cog$trial_type)


write.csv(all_cog,file="ParseANT.csv")
