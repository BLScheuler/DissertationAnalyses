# Exploratory Step 5
# ANT cortisol comparison

### SET UP THE DATA

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)
library(dplyr)

#set working directory

#load capacity coefficient data
ANTdata <- read.csv("ANTparameters.csv")
ANTdata <- ANTdata %>%
  mutate(Week = as.integer(week)) %>%  #fix capitalization and make interger for left join
  mutate(Subject = as.integer(subject)) %>% #fix capitalization and make interger for left join
  select(Subject, Week, sda, rd) #Filter so only keeping columns of interest


#Create subject map for matching participant and subject number
subject_map <- c(
  "Eclipse" = 1,
  "Granite" = 2,
  "Marigold" = 3,
  "Prism" = 4,
  "Quartz" = 5
)

#load self-report data
SRdata <- read.csv("CleanSelfRep.csv")

#rename to match, and only keep data I'm analyzing here
SRdata <- SRdata %>%
  mutate(Subject = subject_map[Participant]) %>%
  mutate(Week = as.integer(Week)) %>%
  filter(Subject %in% 1:5) %>%
  select(Subject, Week, PSS_Total, SF_EmoWB)

#load cortisol data
Cortdata <- read.csv("CortisolResults.csv")
#Marigold missing from Cortdata, so add them in manually (no successful cortisol pulls)
# Create NA rows for Marigold across all 3 weeks
Marigold_missing <- data.frame(
  Participant = "Marigold",
  Week = 1:3,
  EntireSlope = NA,
  DiurnalSlope = NA,
  AUCg = NA
)
Cortdata <- bind_rows(Cortdata, Marigold_missing)

#rename to match, and only keep data I'm analyzing here
Cortdata <- Cortdata %>%
  mutate(Subject = subject_map[Participant]) %>%
  mutate(Week = as.integer(Week)) %>%
  filter(Subject %in% 1:5) %>%
  select(Subject, Week, EntireSlope, DiurnalSlope, AUCg) #Filter so only keeping columns of interest


#Combine in one dataset
AllDat <- Cortdata %>%
  left_join(SRdata, by = c("Subject", "Week")) %>%
  left_join(ANTdata, by = c("Subject", "Week"))

# Make week ordered again, and subject a factor
AllDat <- AllDat %>%
  mutate(Week = factor(Week, levels = 1:3, ordered = TRUE)) %>%
  mutate(Subject = factor(Subject))


#Set up data for modeling
ITbothData <- list(
  part = as.integer(as.factor(AllDat$Subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  DCS = -1*(standardize(AllDat$DiurnalSlope)), #flipped to match direction of stress in model
  AUC = standardize(AllDat$AUCg), #standard normalization
  IT = standardize((AllDat$sda)/(AllDat$rd)) #do Standard Normalization
)


#Set seed for reproducibility
set.seed(123)

### ULAM MODELS

#Both
# PARTRIAL MEDIAION MODEL (TIME -> ATTENTIONAL CONTROL <- STRESS; TIME -> STRESS)
m_BothIT <- ulam(
  alist(
    ## TIME -> ATTENTIONAL CONTROL <- STRESS
    #distribution for IT parameter 
    IT ~ dnorm(mu_IT, sigma_IT),
    #Set up participant change in IT over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_IT <- p_IT[part] + p_ITwk2[part]*Week2 + p_ITwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_IT, p_ITwk2, p_ITwk3)[part] ~ multi_normal( c(a_IT,b1_ITwk, b2_ITwk) , Rho_IT, sigma_indIT),
    
    #IT priors
    a_IT ~ normal(0, 1), 
    b1_ITwk ~ normal(0, 1),
    b2_ITwk ~ normal(0, 1),
    
    
    #Priors for stress
    b_PSS ~ normal(0, 1),
    b_EmoWB ~ normal(0, 1), 
    b_DCS ~ normal(0, 1), 
    b_AUC ~ normal(0, 1), 
    
    sigma_IT ~ exponential(1),
    sigma_indIT ~ exponential(1),
    Rho_IT ~ lkj_corr(2),
    
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
    
  ), data = ITbothData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_BothIT, depth=3)



#Set up data for modeling
ITdcsData <- list(
  part = as.integer(as.factor(AllDat$Subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  DCS = -1*(standardize(AllDat$DiurnalSlope)), #flipped to match direction of stress in model
  IT = standardize((AllDat$sda)/(AllDat$rd)) #do Standard Normalization
)

#DCS
# PARTRIAL MEDIAION MODEL (TIME -> ATTENTIONAL CONTROL <- STRESS; TIME -> STRESS)
m_DCSIT <- ulam(
  alist(
    ## TIME -> ATTENTIONAL CONTROL <- STRESS
    #distribution for IT parameter 
    IT ~ dnorm(mu_IT, sigma_IT),
    #Set up participant change in IT over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_IT <- p_IT[part] + p_ITwk2[part]*Week2 + p_ITwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_IT, p_ITwk2, p_ITwk3)[part] ~ multi_normal( c(a_IT,b1_ITwk, b2_ITwk) , Rho_IT, sigma_indIT),
    
    #IT priors
    a_IT ~ normal(0, 1), 
    b1_ITwk ~ normal(0, 1),
    b2_ITwk ~ normal(0, 1),
    
    
    #Priors for stress
    b_PSS ~ normal(0, 1),
    b_EmoWB ~ normal(0, 1), 
    b_DCS ~ normal(0, 1), 
    
    sigma_IT ~ exponential(1),
    sigma_indIT ~ exponential(1),
    Rho_IT ~ lkj_corr(2),
    
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
    
  ), data = ITdcsData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_DCSIT, depth=3)


#Set up data for modeling
ITaucData <- list(
  part = as.integer(as.factor(AllDat$Subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  AUC = standardize(AllDat$AUCg), #standard normalization
  IT = standardize((AllDat$sda)/(AllDat$rd)) #do Standard Normalization
)

#AUC
# PARTRIAL MEDIAION MODEL (TIME -> ATTENTIONAL CONTROL <- STRESS; TIME -> STRESS)
m_AUCIT <- ulam(
  alist(
    ## TIME -> ATTENTIONAL CONTROL <- STRESS
    #distribution for IT parameter 
    IT ~ dnorm(mu_IT, sigma_IT),
    #Set up participant change in IT over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_IT <- p_IT[part] + p_ITwk2[part]*Week2 + p_ITwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_AUC*AUC,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_IT, p_ITwk2, p_ITwk3)[part] ~ multi_normal( c(a_IT,b1_ITwk, b2_ITwk) , Rho_IT, sigma_indIT),
    
    #IT priors
    a_IT ~ normal(0, 1), 
    b1_ITwk ~ normal(0, 1),
    b2_ITwk ~ normal(0, 1),
    
    
    #Priors for stress
    b_PSS ~ normal(0, 1),
    b_EmoWB ~ normal(0, 1), 
    b_AUC ~ normal(0, 1), 
    
    sigma_IT ~ exponential(1),
    sigma_indIT ~ exponential(1),
    Rho_IT ~ lkj_corr(2),
    
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
    
  ), data = ITaucData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_AUCIT, depth=3)


### MODEL COMPARISON

#Compare WAICs
ITcompare <- compare(m_BothIT, m_DCSIT, m_AUCIT, func = WAIC)
ITcomp_df <- as.data.frame(ITcompare)
ITcomp_df$model <- rownames(ITcomp_df)

#Renaming to look cleaner on plots
ITcomp_df$model[grepl("m_BothIT", ITcomp_df$model)] <- "Both"
ITcomp_df$model[grepl("m_DCSIT", ITcomp_df$model)] <- "DCS"
ITcomp_df$model[grepl("m_AUCIT", ITcomp_df$model)] <- "AUC"

#move model column to front
ITcomp_df <- ITcomp_df[, c("model", setdiff(names(ITcomp_df), "model"))]

# Save to CSV
write.csv(ITcomp_df, "Output/Expl_IT_Cort_ModelComp.csv", row.names = FALSE)
