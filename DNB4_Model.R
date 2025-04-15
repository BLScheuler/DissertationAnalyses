#DNB Step 4
#Modeling the DNB data with the SFT package

#Set working directory

library(SuppDists)
library(sft)
library(dplyr)


ALLdata <- read.csv("dnb.csv")
#number obs 6480

CleanDNB <- read.csv("CleanDNB.csv")
#number obs 6480

#change participant names to subject names, to match other list
ALLdata <- ALLdata %>%
  mutate(subject = recode(subject,
                          "obsidian" = "1",
                          "redwood" = "2",
                          "solstice" = "3",
                          "starlight" = "4",
                          "tundra" = "5"))

#Remove the parse version of rt's, since the CleanDNB is more accurate 
ALLdata <- ALLdata[, -8]


#Need to add trial_type to ALLdata
ALLdata <- ALLdata %>%
  mutate(trial_type = case_when(
    visual != "None" & auditory == "None" ~ "visual",
    visual == "None" & auditory != "None" ~ "audio",
    visual != "None" & auditory != "None" ~ "dual",
    TRUE ~ NA_character_  # just in case something unexpected slips through
  ))


#Add accuracy and rt columns 
mergedData <- merge(ALLdata,
                    CleanDNB[, c("subject", "week", "block_n", "trial_n", "trial_type", "rt", "accuracy")],
                    by.x = c("subject", "week", "block_no", "trial_no", "trial_type"),
                    by.y = c("subject", "week", "block_n", "trial_n", "trial_type"),
                    all.x = TRUE)


#Now need to set up visual and auditory for channels 
#Make it based on AND (with the negatives) to make the target match our design

mergedData <- mergedData %>%
  group_by(subject, week, block_no, trial_type) %>%
  arrange(trial_no, .by_group = TRUE) %>%
  mutate(
    visual_lag2 = lag(visual, 2),
    
    ## Create VisChan
    VisChan = case_when(
      visual == "None" ~ "0",
      visual == visual_lag2 ~ "-2",
      visual != visual_lag2 ~ "2",
      TRUE ~ NA_character_
    ),
    
    ## Create AudChan 
    AudChan = case_when(
      
      # VISUAL ONLY TRIALS
      auditory == "None" ~ "0",
      
      # AUDIO ONLY TRIALS
      trial_type == "audio" & choice == "left" & accuracy == 1 ~ "2",   # correct allow > mismatch
      trial_type == "audio" & choice == "right" & accuracy == 1 ~ "-2", # correct block > match
      trial_type == "audio" & choice == "left" & accuracy == 0 ~ "-2",  # incorrect allow > match
      trial_type == "audio" & choice == "right" & accuracy == 0 ~ "2",  # incorrect block > mismatch
      
      # DUAL TRIALS, where visual ≠ visual_lag2 (so audio is the deciding factor) 
      trial_type == "dual" & visual != visual_lag2 & choice == "left" & accuracy == 1 ~ "2",   # correct allow > audio mismatch
      trial_type == "dual" & visual != visual_lag2 & choice == "right" & accuracy == 1 ~ "-2", # correct block > audio match
      trial_type == "dual" & visual != visual_lag2 & choice == "left" & accuracy == 0 ~ "-2",  # incorrect allow > audio match
      trial_type == "dual" & visual != visual_lag2 & choice == "right" & accuracy == 0 ~ "2",  # incorrect block > audio mismatch
      
      # Problem Children DUAL TRIALS, where visual == visual_lag2 > we don't know what audio is happening
      trial_type == "dual" & visual == visual_lag2 ~ NA_character_,
      
      # fill in other ambiguous cases as NA for safety
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# Subset of data for the allow decision

allow_dat <- subset(mergedData, 
                    (choice=="left" & accuracy == 1 & trial_type=="dual") |
                      (choice=="right" & accuracy == 0 & trial_type =="dual") |
                      (choice=="left" & accuracy == 1 & trial_type=="audio") |
                      (choice=="right" & accuracy == 0 & trial_type =="audio") |
                      (choice=="left" & accuracy == 1 & trial_type=="visual") |
                      (choice=="right" & accuracy == 0 & trial_type =="visual"))
#Remove first 2 trials
allow_dat <- allow_dat[allow_dat$trial_no > 2, ]
allow_dat <- allow_dat[,c("subject", "week", "AudChan", "VisChan", "accuracy", "rt")]
allow_dat <- as.data.frame(allow_dat)
names(allow_dat) <- c("Subject", "Condition", "Channel1", "Channel2", "Correct", "RT") 

allow_dat$Channel1 <- as.numeric(levels(allow_dat$Channel1)[allow_dat$Channel1])
allow_dat$Channel2 <- as.numeric(levels(allow_dat$Channel2)[allow_dat$Channel2])

assess_allow <- assessmentGroup(allow_dat, detection = FALSE, correct = TRUE, fast = TRUE, stopping.rule = "AND")
cap_allow <- capacityGroup(allow_dat, acc.cutoff=.2, stopping.rule="AND")

#Create a dataframe with all the Capacity Coefficients
cap_z_allow <- matrix(NA, nrow=15, ncol=3)
for (sj in 1:5) { 
  for (cn in 1:3) {
    i <- 5 *(cn-1) + sj
    cz <- cap_allow$capacity[[i]]$Ctest$statistic
    cap_z_allow[i,] <- c(sj, cn, cz)
  }
}
cap_z_allow <- as.data.frame(cap_z_allow)
names(cap_z_allow) <- c("Subject", "Week", "Cz")
cap_z_allow$Subject <- factor(cap_z_allow$Subject)
cap_z_allow$Week <- factor(cap_z_allow$Week)


# Save to CSV
write.csv(cap_z_allow, "DNB_CapCoef.csv", row.names = FALSE)


