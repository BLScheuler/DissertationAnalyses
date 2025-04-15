# Exploratory Step 6
# DNB cortisol comparison


### SET UP THE DATA

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)
library(dplyr)

#set working directory

#load capacity coefficient data
DNBdata <- read.csv("DNB_CapCoef.csv")

#Create subject map for matching participant and subject number
subject_map <- c(
  "Obsidian" = 1,
  "Redwood" = 2,
  "Solstice" = 3,
  "Starlight" = 4,
  "Tundra" = 5
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
#Redwood and Solstice are missing from Cortdata, so add them in manually (no successful cortisol pulls)
# Create NA rows for Marigold across all 3 weeks
Redwood_missing <- data.frame(
  Participant = "Redwood",
  Week = 1:3,
  EntireSlope = NA,
  DiurnalSlope = NA,
  AUCg = NA
)
Solstice_missing <- data.frame(
  Participant = "Solstice",
  Week = 1:3,
  EntireSlope = NA,
  DiurnalSlope = NA,
  AUCg = NA
)
Cortdata <- bind_rows(Cortdata, Redwood_missing, Solstice_missing)

#rename to match, and only keep data I'm analyzing here
Cortdata <- Cortdata %>%
  mutate(Subject = subject_map[Participant]) %>%
  mutate(Week = as.integer(Week)) %>%
  filter(Subject %in% 1:5) %>%
  select(Subject, Week, EntireSlope, DiurnalSlope, AUCg) #Filter so only keeping columns of interest


#Combine in one dataset
AllDat <- Cortdata %>%
  left_join(SRdata, by = c("Subject", "Week")) %>%
  left_join(DNBdata, by = c("Subject", "Week"))

# Make week ordered again, and subject a factor
AllDat <- AllDat %>%
  mutate(Week = factor(Week, levels = 1:3, ordered = TRUE)) %>%
  mutate(Subject = factor(Subject))


#Set up data for modeling
WCbothData <- list(
  part = as.integer(as.factor(AllDat$Subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  DCS = -1*(standardize(AllDat$DiurnalSlope)), #flipped to match direction of stress in model
  AUC = (standardize(AllDat$AUCg)), #flipped to match direction of stress in model
  WC = DNBdata$Cz
)


#Set seed for reproducibility
set.seed(123)

### ULAM MODELS

#Both
# PARTRIAL MEDIAION MODEL (TIME -> WORKLOAD CAPACITY <- STRESS; TIME -> STRESS)
m_BothWC <- ulam(
  alist(
    ## TIME -> WORKLOAD CAPACITY <- STRESS
    #distribution for WC parameter 
    WC ~ dnorm(mu_WC, sigma_WC),
    #Set up participant change in WC over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_WC <- p_WC[part] + p_WCwk2[part]*Week2 + p_WCwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_WC, p_WCwk2, p_WCwk3)[part] ~ multi_normal( c(a_WC,b1_WCwk, b2_WCwk) , Rho_WC, sigma_indWC),
    
    #WC priors
    a_WC ~ normal(3, 0.75),
    b1_WCwk ~ normal(0, 0.5),
    b2_WCwk ~ normal(0, 0.5),
    
    #Priors for stress
    b_PSS ~ normal(0,.1),
    b_EmoWB ~ normal(0,.1), 
    b_DCS ~ normal(0, .1), 
    b_AUC ~ normal(0, .1),
    
    sigma_WC ~ exponential(1),
    sigma_indWC ~ exponential(1),
    Rho_WC ~ lkj_corr(2),
    
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
    
  ), data = WCbothData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_BothWC, depth=3)



#Set up data for modeling
WCdcsData <- list(
  part = as.integer(as.factor(AllDat$Subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  DCS = -1*(standardize(AllDat$DiurnalSlope)), #flipped to match direction of stress in model
  WC = DNBdata$Cz
)

#DCS
# PARTRIAL MEDIAION MODEL (TIME -> WORKLOAD CAPACITY <- STRESS; TIME -> STRESS)
m_DCSWC <- ulam(
  alist(
    ## TIME -> WORKLOAD CAPACITY <- STRESS
    #distribution for WC parameter 
    WC ~ dnorm(mu_WC, sigma_WC),
    #Set up participant change in WC over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_WC <- p_WC[part] + p_WCwk2[part]*Week2 + p_WCwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_WC, p_WCwk2, p_WCwk3)[part] ~ multi_normal( c(a_WC,b1_WCwk, b2_WCwk) , Rho_WC, sigma_indWC),
    
    #WC priors
    a_WC ~ normal(3, 0.75),
    b1_WCwk ~ normal(0, 0.5),
    b2_WCwk ~ normal(0, 0.5),
    
    #Priors for stress
    b_PSS ~ normal(0,.1),
    b_EmoWB ~ normal(0,.1), 
    b_DCS ~ normal(0, .1), 
    
    sigma_WC ~ exponential(1),
    sigma_indWC ~ exponential(1),
    Rho_WC ~ lkj_corr(2),
    
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
    
  ), data = WCdcsData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_DCSWC, depth=3)



#Set up data for modeling
WCaucData <- list(
  part = as.integer(as.factor(AllDat$Subject)),
  Week3 = 1*(AllDat$Week == 3),
  Week2 = 1*(AllDat$Week == 2),
  PSS = (AllDat$PSS_Total)/50, #scaled to percent of maximum
  EmoWB = -1*((AllDat$SF_EmoWB)/100), #scaled to percent of maximum; flipped to match direction of stress in model
  AUC = (standardize(AllDat$AUCg)), #flipped to match direction of stress in model
  WC = DNBdata$Cz
)


#AUC
# PARTRIAL MEDIAION MODEL (TIME -> WORKLOAD CAPACITY <- STRESS; TIME -> STRESS)
m_AUCWC <- ulam(
  alist(
    ## TIME -> WORKLOAD CAPACITY <- STRESS
    #distribution for WC parameter 
    WC ~ dnorm(mu_WC, sigma_WC),
    #Set up participant change in WC over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_WC <- p_WC[part] + p_WCwk2[part]*Week2 + p_WCwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_AUC*AUC,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_WC, p_WCwk2, p_WCwk3)[part] ~ multi_normal( c(a_WC,b1_WCwk, b2_WCwk) , Rho_WC, sigma_indWC),
    
    #WC priors
    a_WC ~ normal(3, 0.75),
    b1_WCwk ~ normal(0, 0.5),
    b2_WCwk ~ normal(0, 0.5),
    
    #Priors for stress
    b_PSS ~ normal(0,.1),
    b_EmoWB ~ normal(0,.1), 
    b_AUC ~ normal(0, .1), 
    
    sigma_WC ~ exponential(1),
    sigma_indWC ~ exponential(1),
    Rho_WC ~ lkj_corr(2),
    
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
    
  ), data = WCaucData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_AUCWC, depth=3)


### MODEL COMPARISON

#Compare WAICs
WCcompare <- compare(m_BothWC, m_DCSWC, m_AUCWC, func = WAIC)
WCcomp_df <- as.data.frame(WCcompare)
WCcomp_df$model <- rownames(WCcomp_df)

#Renaming to look cleaner on plots
WCcomp_df$model[grepl("m_BothWC", WCcomp_df$model)] <- "Both"
WCcomp_df$model[grepl("m_DCSWC", WCcomp_df$model)] <- "DCS"
WCcomp_df$model[grepl("m_AUCWC", WCcomp_df$model)] <- "AUC"

#move model column to front
WCcomp_df <- WCcomp_df[, c("model", setdiff(names(WCcomp_df), "model"))]

# Save to CSV
write.csv(WCcomp_df, "Output/Expl_WC_Cort_ModelComp.csv", row.names = FALSE)
