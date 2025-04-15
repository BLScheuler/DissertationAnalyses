# Exploratory Step 8
# RawRT cortisol comparison

### SET UP THE DATA

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)
library(dplyr)

#Set seed for reproducibility
set.seed(123)  

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
ECbothData <- list(
  part = as.integer(as.factor(AllDat$subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  DCS = -1*(standardize(AllDat$DiurnalSlope)), #flipped to match direction of stress in model
  AUC = standardize(AllDat$AUCg), #standard normalization
  RT = AllDat$meanRT
)



### ULAM MODELS

#Both
# PARTRIAL MEDIAION MODEL (TIME -> EXECUTIVE CONTROL <- STRESS; TIME -> STRESS)
m_BothEC <- ulam(
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
    b_PSS ~ normal(0, .1),
    b_EmoWB ~ normal(0, .1), 
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
    
    a_DCS ~ normal(0, 1),  
    b1_DCSwk ~ normal(0, 1),
    b2_DCSwk ~ normal(0, 1),
    
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
    
  ), data = ECbothData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_BothEC, depth=3)




#Set up data for modeling
ECdcsData <- list(
  part = as.integer(as.factor(AllDat$subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  DCS = -1*(standardize(AllDat$DiurnalSlope)), #flipped to match direction of stress in model
  RT = AllDat$meanRT
)

#DCS
# PARTRIAL MEDIAION MODEL (TIME -> EXECUTIVE CONTROL <- STRESS; TIME -> STRESS)
m_DCSEC <- ulam(
  alist(
    ## TIME -> EXECUTIVE CONTROL <- STRESS
    #distribution for EC parameter 
    RT ~ dnorm(mu_RT, sigma_RT),
    #Set up participant change in RT over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_RT <- p_RT[part] + p_RTwk2[part]*Week2 + p_RTwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS,
    
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
    
    sigma_RT ~ exponential(1),
    sigma_indRT ~ exponential(1),
    Rho_RT ~ lkj_corr(2),
    
    ##TIME -> STRESS
    PSS ~ dnorm(mu_PSS, sigma_distPSS),
    EmoWB ~ dnorm(mu_EmoWB, sigma_distEMO),
    DCS ~ dnorm(mu_DCS, sigma_distDCS),
    
    #Set up participant change over time for stress measures
    mu_PSS <- p_PSS[part] + p_PSSwk2[part]*Week2 + p_PSSwk3[part]*Week3,
    mu_EmoWB <- p_EmoWB[part] + p_EMOwk2[part]*Week2 + p_EMOwk3[part]*Week3,
    mu_DCS <- p_DCS[part] + p_DCSwk2[part]*Week2 + p_DCSwk3[part]*Week3,
    
    c(p_PSS,p_PSSwk2, p_PSSwk3)[part] ~ multi_normal( c(a_PSS,b1_PSSwk, b2_PSSwk) , Rho_PSS , sigma_PSS),
    c(p_EmoWB,p_EMOwk2, p_EMOwk3)[part] ~ multi_normal( c(a_EMO,b1_EMOwk, b2_EMOwk) , Rho_EMO , sigma_EMO),
    c(p_DCS,p_DCSwk2, p_DCSwk3)[part] ~ multi_normal( c(a_DCS,b1_DCSwk, b2_DCSwk) , Rho_DCS , sigma_DCS), 
    
    a_PSS ~ normal(.5,.1),
    b1_PSSwk ~ normal(0,.1),
    b2_PSSwk ~ normal(0,.1),
    
    a_EMO ~ normal(.5, .1),
    b1_EMOwk ~ normal(0,.1),
    b2_EMOwk ~ normal(0,.1),
    
    a_DCS ~ normal(0, 1),  
    b1_DCSwk ~ normal(0, 1),
    b2_DCSwk ~ normal(0, 1),
    
    sigma_PSS ~ exponential(1),
    sigma_EMO ~ exponential(1),
    sigma_DCS ~ exponential(1),
    
    Rho_PSS ~ lkj_corr(2),
    Rho_EMO ~ lkj_corr(2),
    Rho_DCS ~ lkj_corr(2),
    
    sigma_distPSS ~ exponential(1),
    sigma_distEMO ~ exponential(1),
    sigma_distDCS ~ exponential(1)
    
  ), data = ECdcsData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_DCSEC, depth=3)



#Set up data for modeling
ECaucData <- list(
  part = as.integer(as.factor(AllDat$subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  AUC = standardize(AllDat$AUCg), #standard normalization
  RT = AllDat$meanRT
)

#AUC
# PARTRIAL MEDIAION MODEL (TIME -> EXECUTIVE CONTROL <- STRESS; TIME -> STRESS)
m_AUCEC <- ulam(
  alist(
    ## TIME -> EXECUTIVE CONTROL <- STRESS
    #distribution for EC parameter 
    RT ~ dnorm(mu_RT, sigma_RT),
    #Set up participant change in RT over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_RT <- p_RT[part] + p_RTwk2[part]*Week2 + p_RTwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_AUC*AUC,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_RT,p_RTwk2, p_RTwk3)[part] ~ multi_normal( c(a_RT,b1_RTwk, b2_RTwk) , Rho_RT , sigma_indRT),
    
    #RT prior 
    a_RT ~ normal(.8,.2),
    b1_RTwk ~ normal(0,.1),
    b2_RTwk ~ normal(0,.1),
    
    #Priors for stress
    b_PSS ~ normal(0,.1),
    b_EmoWB ~ normal(0,.1), 
    b_AUC ~ normal(0, .1), 
    
    sigma_RT ~ exponential(1),
    sigma_indRT ~ exponential(1),
    Rho_RT ~ lkj_corr(2),
    
    ##TIME -> STRESS
    PSS ~ dnorm(mu_PSS, sigma_distPSS),
    EmoWB ~ dnorm(mu_EmoWB, sigma_distEMO),
    AUC ~ dnorm(mu_AUC, sigma_distAUC),
    
    #Set up participant change over time for stress measures
    mu_PSS <- p_PSS[part] + p_PSSwk2[part]*Week2 + p_PSSwk3[part]*Week3,
    mu_EmoWB <- p_EmoWB[part] + p_EMOwk2[part]*Week2 + p_EMOwk3[part]*Week3,
    mu_AUC <- p_AUC[part] + p_AUCwk2[part]*Week2 + p_AUCwk3[part]*Week3,
    
    c(p_PSS,p_PSSwk2, p_PSSwk3)[part] ~ multi_normal( c(a_PSS,b1_PSSwk, b2_PSSwk) , Rho_PSS , sigma_PSS),
    c(p_EmoWB,p_EMOwk2, p_EMOwk3)[part] ~ multi_normal( c(a_EMO,b1_EMOwk, b2_EMOwk) , Rho_EMO , sigma_EMO),
    c(p_AUC,p_AUCwk2, p_AUCwk3)[part] ~ multi_normal( c(a_AUC,b1_AUCwk, b2_AUCwk) , Rho_AUC , sigma_AUC), 
    
    a_PSS ~ normal(.5,.1),
    b1_PSSwk ~ normal(0,.1),
    b2_PSSwk ~ normal(0,.1),
    
    a_EMO ~ normal(.5, .1),
    b1_EMOwk ~ normal(0,.1),
    b2_EMOwk ~ normal(0,.1),
    
    a_AUC ~ normal(0, 1),  
    b1_AUCwk ~ normal(0, 1),
    b2_AUCwk ~ normal(0, 1),
    
    sigma_PSS ~ exponential(1),
    sigma_EMO ~ exponential(1),
    sigma_AUC ~ exponential(1),
    
    Rho_PSS ~ lkj_corr(2),
    Rho_EMO ~ lkj_corr(2),
    Rho_AUC ~ lkj_corr(2),
    
    sigma_distPSS ~ exponential(1),
    sigma_distEMO ~ exponential(1),
    sigma_distAUC ~ exponential(1)
    
  ), data = ECaucData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_AUCEC, depth=3)


### MODEL COMPARISON

#Compare WAICs
ECcompare <- compare(m_BothEC, m_DCSEC, m_AUCEC, func = WAIC)
ECcomp_df <- as.data.frame(ECcompare)
ECcomp_df$model <- rownames(ECcomp_df)

#Renaming to look cleaner on plots
ECcomp_df$model[grepl("m_BothEC", ECcomp_df$model)] <- "Both"
ECcomp_df$model[grepl("m_DCSEC", ECcomp_df$model)] <- "DCS"
ECcomp_df$model[grepl("m_AUCEC", ECcomp_df$model)] <- "AUC"

#move model column to front
ECcomp_df <- ECcomp_df[, c("model", setdiff(names(ECcomp_df), "model"))]

# Save to CSV
write.csv(ECcomp_df, "Output/Expl_ECrt_Cort_ModelComp.csv", row.names = FALSE)
