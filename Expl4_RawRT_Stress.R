# Exploratory Step 4
#Raw RT stress mediation


### SET UP THE DATA

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)
library(dplyr)

#set working directory

#load exgSS modeling data
OSARIdata <- read.csv("osari.csv")

#Keep only RTs that are from unsuccessful stop trials
OSARIdata <- OSARIdata[OSARIdata$ss_presented == 1 & OSARIdata$inhibited == 0, ]

#Get the average rawS for each participant and week
MeanRTdata <- OSARIdata %>%
  group_by(subject, week) %>%
  summarize(meanRT = mean(rt, na.rm = TRUE), .groups = "drop")

# Convert from ms to s
MeanRTdata$meanRT <- MeanRTdata$meanRT / 1000

#rename week to match, and only keep columns I'm interested in
MeanRTdata <- MeanRTdata %>%
  rename(Week = week) %>%  # match the capitalization
  mutate(Week = as.integer(Week)) %>% #temporarily make it integer for left_join
  mutate(subject = as.integer(subject)) #temporarily make it integer for left_join

#load self-report data
SRdata <- read.csv("CleanSelfRep.csv")

#Assign correct subject numbers
SRdata <- SRdata %>%
  mutate(subject = recode(Participant,
                          "Cascade" = "1",
                          "Glacier" = "2",
                          "Harbor" = "3",
                          "Horizon" = "4",
                          "Meadow" = "5")) %>%
  mutate(Week = as.integer(Week)) %>% #temporarily make it integer for left_join
  filter(subject %in% c("1", "2", "3", "4", "5")) %>%
  mutate(subject = as.integer(subject)) %>% #temporarily make it integer for left_join
  select(subject, Week, PSS_Total, SF_EmoWB) #Filter so only keeping columns of interest


#load cortisol data
Cortdata <- read.csv("CortisolResults.csv")
#Glacier is missing from Cortdata, so add them in manually (no succesful cortisol pulls)
# Create NA rows for Glacier across all 3 weeks
glacier_missing <- data.frame(
  Participant = "Glacier",
  Week = 1:3,
  EntireSlope = NA,
  DiurnalSlope = NA,
  AUCg = NA
)
Cortdata <- bind_rows(Cortdata, glacier_missing)

#Assign correct subject numbers
Cortdata <- Cortdata %>%
  mutate(subject = recode(Participant,
                          "Cascade" = "1",
                          "Glacier" = "2",
                          "Harbor" = "3",
                          "Horizon" = "4",
                          "Meadow" = "5")) %>%
  mutate(Week = as.integer(Week)) %>% #temporarily make it integer for left_join
  filter(subject %in% c("1", "2", "3", "4", "5")) %>%
  mutate(subject = as.integer(subject)) %>% #temporarily make it integer for left_join
  select(subject, Week, EntireSlope, DiurnalSlope, AUCg) #Filter so only keeping columns of interest



#Combine in one dataset
AllDat <- Cortdata %>%
  left_join(SRdata, by = c("subject", "Week")) %>%
  left_join(MeanRTdata, by = c("subject", "Week"))

# Make week ordered again
AllDat <- AllDat %>%
  mutate(Week = factor(Week, levels = 1:3, ordered = TRUE))

#Set up data for modeling
ECmodData <- list(
  part = as.integer(as.factor(AllDat$subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  DCS = -1*(AllDat$DiurnalSlope), #flipped to match direction of stress in model
  AUC = standardize(AllDat$AUCg), #standard normalization
  RT = AllDat$meanRT
)


### ULAM MODELS

# PARTRIAL MEDIAION MODEL (TIME -> EXECUTIVE CONTROL <- STRESS; TIME -> STRESS)
m_PartMedEC <- ulam(
  alist(
    ## TIME -> EXECUTIVE CONTROL <- STRESS
    #distribution for EC parameter 
    RT ~ dnorm(mu_RT, sigma_RT),
    #Set up participant change in RT over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_RT <- p_RT[part] + p_RTwk2[part]*Week2 + p_RTwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_RT,p_RTwk2, p_RTwk3)[part] ~ multi_normal( c(a_RT,b1_RTwk, b2_RTwk) , Rho_RT , sigma_indRT),
    
    #RT priors 
    a_RT ~ normal(.8,.2),
    b1_RTwk ~ normal(0,.1),
    b2_RTwk ~ normal(0,.1),
    
    #Priors for stress
    b_PSS ~ normal(0,.1),
    b_EmoWB ~ normal(0,.1), 
    b_DCS ~ normal(0, .1), 
    b_AUC ~ normal(0, .1),
    
    sigma_RT ~ exponential(1),
    sigma_indRT ~ exponential(1),
    Rho_RT ~ lkj_corr(2),
    
    ##TIME -> STRESS
    PSS ~ dnorm(mu_PSS, sigma_distPSS),
    EmoWB ~ dnorm(mu_EmoWB, sigma_distEMO),
    DCS ~ dnorm(mu_DCS, sigma_distDCS),
    AUC ~ dnorm(mu_AUC, sigma_distAUC),
    
    #Set up participant change over time for stress measures
    mu_PSS <- p_PSS[part] + p_PSSwk2[part]*Week2 + p_PSSwk3[part]*Week3,
    mu_EmoWB <- p_EmoWB[part] + p_EMOwk2[part]*Week2 + p_EMOwk3[part]*Week3,
    mu_DCS <- p_DCS[part] + p_DCSwk2[part]*Week2 + p_DCSwk3[part]*Week3,
    mu_AUC <- p_AUC[part] + p_AUCwk2[part]*Week2 + p_AUCwk3[part]*Week3,
    
    c(p_PSS,p_PSSwk2, p_PSSwk3)[part] ~ multi_normal( c(a_PSS,b1_PSSwk, b2_PSSwk) , Rho_PSS , sigma_PSS),
    c(p_EmoWB,p_EMOwk2, p_EMOwk3)[part] ~ multi_normal( c(a_EMO,b1_EMOwk, b2_EMOwk) , Rho_EMO , sigma_EMO),
    c(p_DCS,p_DCSwk2, p_DCSwk3)[part] ~ multi_normal( c(a_DCS,b1_DCSwk, b2_DCSwk) , Rho_DCS , sigma_DCS), 
    c(p_AUC,p_AUCwk2, p_AUCwk3)[part] ~ multi_normal( c(a_AUC,b1_AUCwk, b2_AUCwk) , Rho_AUC , sigma_AUC),
    
    a_PSS ~ normal(.5,.1),
    b1_PSSwk ~ normal(0,.1),
    b2_PSSwk ~ normal(0,.1),
    
    a_EMO ~ normal(.5, .1),
    b1_EMOwk ~ normal(0,.1),
    b2_EMOwk ~ normal(0,.1),
    
    a_DCS ~ normal(0.5, .1),  
    b1_DCSwk ~ normal(0,.1),
    b2_DCSwk ~ normal(0,.1),
    
    a_AUC ~ normal(0, 1),  
    b1_AUCwk ~ normal(0, 1),
    b2_AUCwk ~ normal(0, 1),
    
    sigma_PSS ~ exponential(1),
    sigma_EMO ~ exponential(1),
    sigma_DCS ~ exponential(1),
    sigma_AUC ~ exponential(1),
    
    Rho_PSS ~ lkj_corr(2),
    Rho_EMO ~ lkj_corr(2),
    Rho_DCS ~ lkj_corr(2),
    Rho_AUC ~ lkj_corr(2),
    
    sigma_distPSS ~ exponential(1),
    sigma_distEMO ~ exponential(1),
    sigma_distDCS ~ exponential(1),
    sigma_distAUC ~ exponential(1)
    
  ), data = ECmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_PartMedEC, depth=3)




#FULL MEDIATION MODEL (TIME -> STRESS -> EXECUTIVE CONTROL)
m_FullMedEC <- ulam(
  alist(
    ## STRESS -> EXECUTIVE CONTROL
    #distribution for EC parameter 
    RT ~ dnorm(mu_RT, sigma_RT),
    #Set up participant change over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3
    mu_RT <- p_RT[part] + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    
    #RT priors
    p_RT[part] ~ normal(a_RT, sigma_indRT),
    a_RT ~ normal(.8,.2),
    
    #Priors for stress
    b_PSS ~ normal(0,.1),
    b_EmoWB ~ normal(0,.1), 
    b_DCS ~ normal(0, .1), 
    b_AUC ~ normal(0, .1),
    
    sigma_RT ~ exponential(1),
    sigma_indRT ~ exponential(1),
    
    ## TIME -> STRESS
    PSS ~ dnorm(mu_PSS, sigma_distPSS),
    EmoWB ~ dnorm(mu_EmoWB, sigma_distEMO),
    DCS ~ dnorm(mu_DCS, sigma_distDCS),
    AUC ~ dnorm(mu_AUC, sigma_distAUC),
    
    #Set up participant change over time for stress measures
    mu_PSS <- p_PSS[part] + p_PSSwk2[part]*Week2 + p_PSSwk3[part]*Week3,
    mu_EmoWB <- p_EmoWB[part] + p_EMOwk2[part]*Week2 + p_EMOwk3[part]*Week3,
    mu_DCS <- p_DCS[part] + p_DCSwk2[part]*Week2 + p_DCSwk3[part]*Week3,
    mu_AUC <- p_AUC[part] + p_AUCwk2[part]*Week2 + p_AUCwk3[part]*Week3,
    
    c(p_PSS,p_PSSwk2, p_PSSwk3)[part] ~ multi_normal( c(a_PSS,b1_PSSwk, b2_PSSwk) , Rho_PSS , sigma_PSS),
    c(p_EmoWB,p_EMOwk2, p_EMOwk3)[part] ~ multi_normal( c(a_EMO,b1_EMOwk, b2_EMOwk) , Rho_EMO , sigma_EMO),
    c(p_DCS,p_DCSwk2, p_DCSwk3)[part] ~ multi_normal( c(a_DCS,b1_DCSwk, b2_DCSwk) , Rho_DCS , sigma_DCS), 
    c(p_AUC,p_AUCwk2, p_AUCwk3)[part] ~ multi_normal( c(a_AUC,b1_AUCwk, b2_AUCwk) , Rho_AUC , sigma_AUC), 
    
    a_PSS ~ normal(.5,.1),
    b1_PSSwk ~ normal(0,.1),
    b2_PSSwk ~ normal(0,.1),
    
    a_EMO ~ normal(.5, .1),
    b1_EMOwk ~ normal(0,.1),
    b2_EMOwk ~ normal(0,.1),
    
    a_DCS ~ normal(0.5, .1),  
    b1_DCSwk ~ normal(0,.1),
    b2_DCSwk ~ normal(0,.1),
    
    a_AUC ~ normal(0, 1),  
    b1_AUCwk ~ normal(0, 1),
    b2_AUCwk ~ normal(0, 1),
    
    sigma_PSS ~ exponential(1),
    sigma_EMO ~ exponential(1),
    sigma_DCS ~ exponential(1),
    sigma_AUC ~ exponential(1),
    
    Rho_PSS ~ lkj_corr(2),
    Rho_EMO ~ lkj_corr(2),
    Rho_DCS ~ lkj_corr(2),
    Rho_AUC ~ lkj_corr(2),
    
    sigma_distPSS ~ exponential(1),
    sigma_distEMO ~ exponential(1),
    sigma_distDCS ~ exponential(1),
    sigma_distAUC ~ exponential(1)
    
  ), data = ECmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_FullMedEC, depth=3)




#NO MEDIATION MODEL (TIME -> EXECUTIVE CONTROL; TIME -> STRESS)
m_TimeEC <- ulam(
  alist(
    ## TIME -> EXECUTIVE CONTROL
    #distribution for EC parameter 
    RT ~ dnorm(mu_RT, sigma_RT),
    #Set up participant change over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3
    mu_RT <- p_RT[part] + p_RTwk2[part]*Week2 + p_RTwk3[part]*Week3,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_RT,p_RTwk2, p_RTwk3)[part] ~ multi_normal( c(a_RT,b1_RTwk, b2_RTwk) , Rho_RT , sigma_indRT),
    
    #RT prior 
    a_RT ~ normal(.8,.2),
    b1_RTwk ~ normal(0,.1),
    b2_RTwk ~ normal(0,.1),
    
    
    sigma_RT ~ exponential(1),
    sigma_indRT ~ exponential(1),
    Rho_RT ~ lkj_corr(2),
    
    ##TIME -> STRESS
    PSS ~ dnorm(mu_PSS, sigma_distPSS),
    EmoWB ~ dnorm(mu_EmoWB, sigma_distEMO),
    DCS ~ dnorm(mu_DCS, sigma_distDCS),
    AUC ~ dnorm(mu_AUC, sigma_distAUC),
    
    #Set up participant change over time for stress measures
    mu_PSS <- p_PSS[part] + p_PSSwk2[part]*Week2 + p_PSSwk3[part]*Week3,
    mu_EmoWB <- p_EmoWB[part] + p_EMOwk2[part]*Week2 + p_EMOwk3[part]*Week3,
    mu_DCS <- p_DCS[part] + p_DCSwk2[part]*Week2 + p_DCSwk3[part]*Week3,
    mu_AUC <- p_AUC[part] + p_AUCwk2[part]*Week2 + p_AUCwk3[part]*Week3,
    
    c(p_PSS,p_PSSwk2, p_PSSwk3)[part] ~ multi_normal( c(a_PSS,b1_PSSwk, b2_PSSwk) , Rho_PSS , sigma_PSS),
    c(p_EmoWB,p_EMOwk2, p_EMOwk3)[part] ~ multi_normal( c(a_EMO,b1_EMOwk, b2_EMOwk) , Rho_EMO , sigma_EMO),
    c(p_DCS,p_DCSwk2, p_DCSwk3)[part] ~ multi_normal( c(a_DCS,b1_DCSwk, b2_DCSwk) , Rho_DCS , sigma_DCS), 
    c(p_AUC,p_AUCwk2, p_AUCwk3)[part] ~ multi_normal( c(a_AUC,b1_AUCwk, b2_AUCwk) , Rho_AUC , sigma_AUC),
    
    a_PSS ~ normal(.5,.1),
    b1_PSSwk ~ normal(0,.1),
    b2_PSSwk ~ normal(0,.1),
    
    a_EMO ~ normal(.5, .1),
    b1_EMOwk ~ normal(0,.1),
    b2_EMOwk ~ normal(0,.1),
    
    a_DCS ~ normal(0.5, .1),  
    b1_DCSwk ~ normal(0,.1),
    b2_DCSwk ~ normal(0,.1),
    
    a_AUC ~ normal(0, 1),  
    b1_AUCwk ~ normal(0, 1),
    b2_AUCwk ~ normal(0, 1),
    
    sigma_PSS ~ exponential(1),
    sigma_EMO ~ exponential(1),
    sigma_DCS ~ exponential(1),
    sigma_AUC ~ exponential(1),
    
    Rho_PSS ~ lkj_corr(2),
    Rho_EMO ~ lkj_corr(2),
    Rho_DCS ~ lkj_corr(2),
    Rho_AUC ~ lkj_corr(2),
    
    sigma_distPSS ~ exponential(1),
    sigma_distEMO ~ exponential(1),
    sigma_distDCS ~ exponential(1),
    sigma_distAUC ~ exponential(1)
    
  ), data = ECmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_TimeEC, depth=3)


### POSTERIOR PREDICTIVE PLOTS

#How well did the Partial Mediation Model approximate the posterior distribution?

postP <- extract.samples(m_PartMedEC)

# Number of samples and observations
n_samples <- length(postP$a_RT)
n_obs <- length(ECmodData$RT)

# Storage for predicted means
muP_RT_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ECmodData$part[i]
  w2 <- ECmodData$Week2[i]
  w3 <- ECmodData$Week3[i]
  
  # Pull stress predictors (may be NA)
  pss <- ECmodData$PSS[i]
  emo <- ECmodData$EmoWB[i]
  dcs <- ECmodData$DCS[i]
  auc <- ECmodData$AUC[i]
  
  for (s in 1:n_samples) {
    base <- postP$p_RT[s, pid]
    wk2 <- postP$p_RTwk2[s, pid] * w2
    wk3 <- postP$p_RTwk3[s, pid] * w3
    
    stress <- 0
    if (!is.na(pss)) stress <- stress + postP$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postP$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postP$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postP$b_AUC[s] * auc
    
    muP_RT_manual[s, i] <- base + wk2 + wk3 + stress
  }
}

muP_mean <- apply(muP_RT_manual, 2, mean)
muP_PI <- apply(muP_RT_manual, 2, PI)


#How well did the Full Mediation Model approximate the posterior distribution?
postF <- extract.samples(m_FullMedEC)

muF_RT_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ECmodData$part[i]
  pss <- ECmodData$PSS[i]
  emo <- ECmodData$EmoWB[i]
  dcs <- ECmodData$DCS[i]
  auc <- ECmodData$AUC[i]
  
  for (s in 1:n_samples) {
    base <- postF$p_RT[s, pid]
    stress <- 0
    if (!is.na(pss)) stress <- stress + postF$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postF$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postF$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postF$b_AUC[s] * auc
    
    muF_RT_manual[s, i] <- base + stress
  }
}

muF_mean <- apply(muF_RT_manual, 2, mean)
muF_PI <- apply(muF_RT_manual, 2, PI)

#How well did the No Mediation Model approximate the posterior distribution?
postN <- extract.samples(m_TimeEC)

muN_RT_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ECmodData$part[i]
  w2 <- ECmodData$Week2[i]
  w3 <- ECmodData$Week3[i]
  
  for (s in 1:n_samples) {
    base <- postN$p_RT[s, pid]
    wk2 <- postN$p_RTwk2[s, pid] * w2
    wk3 <- postN$p_RTwk3[s, pid] * w3
    
    muN_RT_manual[s, i] <- base + wk2 + wk3
  }
}

muN_mean <- apply(muN_RT_manual, 2, mean)
muN_PI <- apply(muN_RT_manual, 2, PI)



## All the plots in one image
# Open PNG device with higher resolution
png("Images/Expl_RT_Combined_PostPredPlot.png", width = 1800, height = 600, res = 150)

# Set up 1 row, 3 columns layout and increase overall text size
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2), cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.6)

# --- Partial Mediation Model ---
plot(muP_mean ~ ECmodData$RT, col = rangi2, ylim = range(muP_PI),
     main = "Partial Mediation Model", xlab = "Observed RT", ylab = "Predicted RT")
abline(a = 0, b = 1, lty = 2)
for (i in 1:length(ECmodData$RT)) {
  lines(rep(ECmodData$RT[i], 2), muP_PI[, i], col = rangi2)
}

# --- Full Mediation Model ---
plot(muF_mean ~ ECmodData$RT, col = rangi2, ylim = range(muF_PI),
     main = "Full Mediation Model", xlab = "Observed RT", ylab = "Predicted RT")
abline(a = 0, b = 1, lty = 2)
for (i in 1:length(ECmodData$RT)) {
  lines(rep(ECmodData$RT[i], 2), muF_PI[, i], col = rangi2)
}

# --- No Mediation Model ---
plot(muN_mean ~ ECmodData$RT, col = rangi2, ylim = range(muN_PI),
     main = "No Mediation Model", xlab = "Observed RT", ylab = "Predicted RT")
abline(a = 0, b = 1, lty = 2)
for (i in 1:length(ECmodData$RT)) {
  lines(rep(ECmodData$RT[i], 2), muN_PI[, i], col = rangi2)
}

# Close the device
dev.off()


### MODEL COMPARISON

#Compare WAICs
RTcompare <- compare(m_PartMedEC, m_FullMedEC, m_TimeEC, func = WAIC)
RTcomp_df <- as.data.frame(RTcompare)
RTcomp_df$model <- rownames(RTcomp_df)

#Renaming to look cleaner on plots
RTcomp_df$model[grepl("m_TimeEC", RTcomp_df$model)] <- "No Mediation"
RTcomp_df$model[grepl("m_PartMedEC", RTcomp_df$model)] <- "Partial Mediation"
RTcomp_df$model[grepl("m_FullMedEC", RTcomp_df$model)] <- "Full Mediation"

#move model column to front
RTcomp_df <- RTcomp_df[, c("model", setdiff(names(RTcomp_df), "model"))]

# Save to CSV
write.csv(RTcomp_df, "Output/Expl_RT_ModelComp.csv", row.names = FALSE)

#Rethinking Plot 
plot(compare(m_PartMedEC, m_FullMedEC, m_TimeEC, func = WAIC))  
#See Statistical Rethinking, Chapter 7, code line 7.29

