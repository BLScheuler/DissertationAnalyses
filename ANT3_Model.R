#ANT Step 3
#Fitting ANT with SSP Model

#set working directory

#Load packages
library(devtools)
library(flankr)
library(coda)
library(dplyr)


#Load the data
ANTdata <- read.csv("CleanANT.csv")

#What are the participant numbers?
subs <- unique(ANTdata$subject)

#Adding in weeks
weeks <- unique(ANTdata$week)

#What starting parameters should be used? 
#using parameters listed on rdrr.io for fitSSP
#must be in order: A, ter, p, rd, sda
parameters = c(0.05, 0.3, 0.4, 0.05, 1.8)

#generate empty matrix to store best-fitting parapmeters in for all participants
#subjectfits <- matrix(0, nrow = length(subs), ncol = length(parameters))
subjectfits <- matrix(0, nrow = 15, ncol = 7)
colnames(subjectfits) <- c("subject", "week", "A", "Ter", "p", "rd", "sda")

#reset looping variable
i = 1

##start the subject fits here
for(currSubject in subs){
  
  #loop for each week (Bry edition)
  for(currWeek in weeks){
    
    #get the current subject's data
    subData <- subset(ANTdata, ANTdata$subject == currSubject & ANTdata$week ==currWeek)
    
    #fit the model
    fit <- fitSSP_fixed(subData, parms = parameters, nTrials = 10000,
                        multipleSubjects = FALSE, fixed = c(FALSE, FALSE, FALSE, FALSE, TRUE))
    
    #now store the best-fitting parameters
    subjectfits[i, ] <- c(currSubject, currWeek, fit$bestParameters) 
    
    #update looping variable
    i <- i + 1
  }}


write.csv(subjectfits, file = "ANTparameters.csv", row.names = FALSE)

