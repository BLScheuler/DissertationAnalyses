#DNB Step 2
#Data cleaning


library(rethinking)
library(plyr)
library(dplyr)


#Set working directory

data <- read.csv("ParseDNB.csv")
#number obs 8010

#Subset for Exp trials only
cleanData <- data %>%
  filter(grepl("Exp", trial_type))
#number obs 8010 > 6480


cleanData$participant <- factor(cleanData$participant, ordered = TRUE)
cleanData$subject <- as.integer(cleanData$participant)
cleanData$week <- factor(cleanData$week)
cleanData$accuracy <- as.numeric(as.logical(cleanData$correct))
cleanData$rt <-cleanData$responseTime
cleanData$response <- factor(cleanData$response)
cleanData$block_n <- factor(cleanData$block_n, ordered=TRUE)
cleanData$trial_n <- factor(cleanData$trial_n, ordered=TRUE)



#Remove first two trials
#Don't do this yet, cause it messes with modeling code
#cleanData <- cleanData[cleanData$trial_n > 2, ]
#number obs 6480 > 6000

#Remove too fast
#cleanData <- cleanData[cleanData$rt >= 0.25, ]
#number obs 6000 to 5920

#Remove non-responses
#Note to self: is.na handles "NA" while != "" handles nothing in the cell (I think)
#cleanData <- cleanData[!is.na(cleanData$rt) & cleanData$rt != "", ]
#number obs stayed 5920



#Rename trial type
cleanData$trial_type[grepl("Audio", cleanData$trial_type)] <- "audio"
cleanData$trial_type[grepl("Visual", cleanData$trial_type)] <- "visual"
cleanData$trial_type[grepl("Dual", cleanData$trial_type)] <- "dual"


#Create CSV
CleanDNB <- cleanData[, c("participant", "subject", "week", "response", "correct", "accuracy", "rt", "trial_type", "trial_n", "block_n")]

# Write to CSV
write.csv(CleanDNB, "CleanDNB.csv", row.names = FALSE)


#Make sure no more than 25% was removed per participant
#Original number of trials per participant = 400 (16 x 25)
#partByweek <- count(cleanData, c("participant", "week"))
#partByweek$freq <- partByweek$freq/400

