#OSARI Step 2
#clean data

library(dplyr)
library(tidyr)

#set working directory

OSARIdata <- read.csv("osari.csv")
#4350 obs

# Apply the same cleaning used during modeling
OSARIdata$subject <- as.numeric(factor(OSARIdata$subject))
OSARIdata <- subset(OSARIdata, rt < 2000)
#4350 obs > 4347 
OSARIdata$rt[OSARIdata$rt >= 1000] <- NA 
OSARIdata$rt <- OSARIdata$rt / 1000  # convert from ms to seconds
OSARIdata$ssd <- as.numeric(OSARIdata$ssd)  # Ensure numeric

#Change the trial type column for clarity
OSARIdata <- OSARIdata %>%
  rename(trial_type = ss_presented) %>%
  mutate(trial_type = ifelse(trial_type == 0, "GO", "STOP")) 



#Drop fast response times
OSARIdata <- OSARIdata %>%
  group_by(subject) %>%
  filter(!is.na(rt) & abs(rt - mean(rt, na.rm = TRUE)) <= 3 * sd(rt, na.rm = TRUE)) %>%
  ungroup()

#4347 obs > 3551

#Create csv
write.csv(OSARIdata, "CleanOSARI.csv", row.names = FALSE)
