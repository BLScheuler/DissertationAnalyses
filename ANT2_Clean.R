#ANT Step 2
#Clean


library(rethinking)
library(plyr)

#Set working directory
#data <- read.csv("ParseANT.csv")


##PREPARE VARIABLES

#Turn participant keywords into integers for later modeling
data$partFact <- as.factor(data$participant)
data$subject <- as.integer(data$partFact)

#Set up rt based on flankr format
data$rt <- data$resp.rt

#Set up accuracy based on flankr format
data$accuracy <- data$resp.corr

#Separating target based on flanker condition
data$congruency[grepl("cong", data$tar)] <- "congruent"
data$congruency[grepl("incong", data$tar)] <- "incongruent"
data$congruency[grepl("neutral", data$tar)] <- "neutral"

#Separating target based on cue condition
data$condition[grepl("both", data$cue)] <- "dual"
data$condition[grepl("centre", data$cue)] <- "center"
data$condition[grepl("upper", data$cue)] <- "spatial"
data$condition[grepl("lower", data$cue)] <- "spatial"
data$condition[grepl("blank", data$cue)] <- "none"


##IF I want to look at differences later...
#Separating target based on arrow direction
data$arrowDir[grepl("Right", data$tar)] <- "Right"
data$arrowDir[grepl("Left", data$tar)] <- "Left"

#Turn arrow direction into integers for later modeling
data$arrowFact <- as.factor(data$arrowDir)
data$arrowInt <- as.integer(data$arrowFact)


#Specifying screen condition
data$screenLoc[grepl("0", data$targOrientation)] <- "Bottom"
data$screenLoc[grepl("180", data$targOrientation)] <- "Top"
#Need help to reconfirm this direction piece in psychoPy code

#Turn screen condition into integers for later modeling (if needed)
data$screenFact <- as.factor(data$screenLoc)
data$screenInt <- as.integer(data$screenFact)


##CLEAN RTs

#Make subset for cleaned responses only
cleanData <- data[data$trial_type == "main", ]
#number obs from 4680 > 4320

#Remove non-responses
#Note to self: is.na handles "NA" while != "" handles nothing in the cell (I think)
cleanData <- cleanData[!is.na(cleanData$resp.rt) & cleanData$resp.rt != "", ]
#number obs from 4320 > 4313

#Remove too fast
cleanData <- cleanData[cleanData$resp.rt >= 0.25, ]
#number obs stayed at 4313

#Remove incorrect (if doing counts)
#cleanData <- cleanData[cleanData$resp.corr == 1, ]
#number obs from 4313 > 4134

#Make sure no more than 25% was removed per participant
#Original number of trials per participant = 288 (3 x 96)
#partByweek <- count(cleanData, c("participant", "week"))
#partByweek$freq <- partByweek$freq/288


#Create CSV
CleanANT <- cleanData[, c("participant", "week", "subject", "condition", "congruency", "accuracy", "rt")]

# Write to CSV
write.csv(CleanANT, "CleanANT.csv", row.names = FALSE)







