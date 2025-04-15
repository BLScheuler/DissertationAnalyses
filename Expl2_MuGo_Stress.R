# Exploratory Step 2
#MuGo stress mediation



### SET UP THE DATA

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)
library(dplyr)

#set working directory

#load exgSS modeling data
load("BeestsSubjMeans.RData")

#rename week to match, and only keep columns I'm interested in
MuGodata <- beests_subj_means %>%
  rename(Week = week) %>%  # match the capitalization
  mutate(Week = as.integer(Week)) %>% #temporarily make it integer for left_join
  select(subject, Week, mu.true)

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
  select(subject, Week, EntireSlope, DiurnalSlope, AUCg) #Filter so only keeping columns of interest





#Combine in one dataset
AllDat <- Cortdata %>%
  left_join(SRdata, by = c("subject", "Week")) %>%
  left_join(MuGodata, by = c("subject", "Week"))

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
  MuGo = AllDat$mu.true
)

#Set seed for reproducibility
set.seed(123)  

### PARTIAL MEDIATION MODEL (TIME -> EXECUTIVE CONTROL <- STRESS; TIME -> STRESS)
m_PartMedEC <- ulam(
  alist(
    ## TIME -> EXECUTIVE CONTROL <- STRESS
    MuGo ~ dnorm(mu_MuGo, sigma_MuGo),
    mu_MuGo <- p_MuGo[part] + p_MuGowk2[part]*Week2 + p_MuGowk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    c(p_MuGo, p_MuGowk2, p_MuGowk3)[part] ~ multi_normal( c(a_MuGo, b1_MuGowk, b2_MuGowk) , Rho_MuGo , sigma_indMuGo),
    
    a_MuGo ~ normal(.824,.04),
    b1_MuGowk ~ normal(0, .1),
    b2_MuGowk ~ normal(0, .1),
    
    #Priors for stress
    b_PSS ~ normal(0, .1),
    b_EmoWB ~ normal(0, .1),
    b_DCS ~ normal(0, .1),
    b_AUC ~ normal(0, .1),
    
    sigma_MuGo ~ exponential(1),
    sigma_indMuGo ~ exponential(1),
    Rho_MuGo ~ lkj_corr(2),
    
    ## TIME -> STRESS
    PSS ~ dnorm(mu_PSS, sigma_distPSS),
    EmoWB ~ dnorm(mu_EmoWB, sigma_distEMO),
    DCS ~ dnorm(mu_DCS, sigma_distDCS),
    AUC ~ dnorm(mu_AUC, sigma_distAUC),
    
    mu_PSS <- p_PSS[part] + p_PSSwk2[part]*Week2 + p_PSSwk3[part]*Week3,
    mu_EmoWB <- p_EmoWB[part] + p_EMOwk2[part]*Week2 + p_EMOwk3[part]*Week3,
    mu_DCS <- p_DCS[part] + p_DCSwk2[part]*Week2 + p_DCSwk3[part]*Week3,
    mu_AUC <- p_AUC[part] + p_AUCwk2[part]*Week2 + p_AUCwk3[part]*Week3,
    
    c(p_PSS, p_PSSwk2, p_PSSwk3)[part] ~ multi_normal( c(a_PSS, b1_PSSwk, b2_PSSwk), Rho_PSS, sigma_PSS),
    c(p_EmoWB, p_EMOwk2, p_EMOwk3)[part] ~ multi_normal( c(a_EMO, b1_EMOwk, b2_EMOwk), Rho_EMO, sigma_EMO),
    c(p_DCS, p_DCSwk2, p_DCSwk3)[part] ~ multi_normal( c(a_DCS, b1_DCSwk, b2_DCSwk), Rho_DCS, sigma_DCS),
    c(p_AUC,p_AUCwk2, p_AUCwk3)[part] ~ multi_normal( c(a_AUC,b1_AUCwk, b2_AUCwk) , Rho_AUC , sigma_AUC),
    
    a_PSS ~ normal(.5, .1),
    b1_PSSwk ~ normal(0, .1),
    b2_PSSwk ~ normal(0, .1),
    
    a_EMO ~ normal(.5, .1),
    b1_EMOwk ~ normal(0, .1),
    b2_EMOwk ~ normal(0, .1),
    
    a_DCS ~ normal(0.5, .1),
    b1_DCSwk ~ normal(0, .1),
    b2_DCSwk ~ normal(0, .1),
    
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


### FULL MEDIATION MODEL (TIME -> STRESS -> EXECUTIVE CONTROL)
m_FullMedEC <- ulam(
  alist(
    ## STRESS -> EXECUTIVE CONTROL
    MuGo ~ dnorm(mu_MuGo, sigma_MuGo),
    mu_MuGo <- p_MuGo[part] + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    p_MuGo[part] ~ normal(a_MuGo, sigma_indMuGo),
    a_MuGo ~ normal(.824, .04),
    
    b_PSS ~ normal(0, .1),
    b_EmoWB ~ normal(0, .1),
    b_DCS ~ normal(0, .1),
    b_AUC ~ normal(0, .1),
    
    sigma_MuGo ~ exponential(1),
    sigma_indMuGo ~ exponential(1),
    
    ## TIME -> STRESS
    PSS ~ dnorm(mu_PSS, sigma_distPSS),
    EmoWB ~ dnorm(mu_EmoWB, sigma_distEMO),
    DCS ~ dnorm(mu_DCS, sigma_distDCS),
    AUC ~ dnorm(mu_AUC, sigma_distAUC),
    
    mu_PSS <- p_PSS[part] + p_PSSwk2[part]*Week2 + p_PSSwk3[part]*Week3,
    mu_EmoWB <- p_EmoWB[part] + p_EMOwk2[part]*Week2 + p_EMOwk3[part]*Week3,
    mu_DCS <- p_DCS[part] + p_DCSwk2[part]*Week2 + p_DCSwk3[part]*Week3,
    mu_AUC <- p_AUC[part] + p_AUCwk2[part]*Week2 + p_AUCwk3[part]*Week3,
    
    c(p_PSS, p_PSSwk2, p_PSSwk3)[part] ~ multi_normal( c(a_PSS, b1_PSSwk, b2_PSSwk), Rho_PSS, sigma_PSS),
    c(p_EmoWB, p_EMOwk2, p_EMOwk3)[part] ~ multi_normal( c(a_EMO, b1_EMOwk, b2_EMOwk), Rho_EMO, sigma_EMO),
    c(p_DCS, p_DCSwk2, p_DCSwk3)[part] ~ multi_normal( c(a_DCS, b1_DCSwk, b2_DCSwk), Rho_DCS, sigma_DCS),
    c(p_AUC,p_AUCwk2, p_AUCwk3)[part] ~ multi_normal( c(a_AUC,b1_AUCwk, b2_AUCwk) , Rho_AUC , sigma_AUC), 
    
    a_PSS ~ normal(.5, .1),
    b1_PSSwk ~ normal(0, .1),
    b2_PSSwk ~ normal(0, .1),
    
    a_EMO ~ normal(.5, .1),
    b1_EMOwk ~ normal(0, .1),
    b2_EMOwk ~ normal(0, .1),
    
    a_DCS ~ normal(0.5, .1),
    b1_DCSwk ~ normal(0, .1),
    b2_DCSwk ~ normal(0, .1),
    
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



### NO MEDIATION MODEL (TIME -> EXECUTIVE CONTROL; TIME -> STRESS)
m_TimeEC <- ulam(
  alist(
    ## TIME -> EXECUTIVE CONTROL
    MuGo ~ dnorm(mu_MuGo, sigma_MuGo),
    mu_MuGo <- p_MuGo[part] + p_MuGowk2[part]*Week2 + p_MuGowk3[part]*Week3,
    
    c(p_MuGo, p_MuGowk2, p_MuGowk3)[part] ~ multi_normal( c(a_MuGo, b1_MuGowk, b2_MuGowk), Rho_MuGo, sigma_indMuGo),
    
    a_MuGo ~ normal(.824, .04),
    b1_MuGowk ~ normal(0, .1),
    b2_MuGowk ~ normal(0, .1),
    
    sigma_MuGo ~ exponential(1),
    sigma_indMuGo ~ exponential(1),
    Rho_MuGo ~ lkj_corr(2),
    
    ## TIME -> STRESS
    PSS ~ dnorm(mu_PSS, sigma_distPSS),
    EmoWB ~ dnorm(mu_EmoWB, sigma_distEMO),
    DCS ~ dnorm(mu_DCS, sigma_distDCS),
    AUC ~ dnorm(mu_AUC, sigma_distAUC),
    
    mu_PSS <- p_PSS[part] + p_PSSwk2[part]*Week2 + p_PSSwk3[part]*Week3,
    mu_EmoWB <- p_EmoWB[part] + p_EMOwk2[part]*Week2 + p_EMOwk3[part]*Week3,
    mu_DCS <- p_DCS[part] + p_DCSwk2[part]*Week2 + p_DCSwk3[part]*Week3,
    mu_AUC <- p_AUC[part] + p_AUCwk2[part]*Week2 + p_AUCwk3[part]*Week3,
    
    c(p_PSS, p_PSSwk2, p_PSSwk3)[part] ~ multi_normal( c(a_PSS, b1_PSSwk, b2_PSSwk), Rho_PSS, sigma_PSS),
    c(p_EmoWB, p_EMOwk2, p_EMOwk3)[part] ~ multi_normal( c(a_EMO, b1_EMOwk, b2_EMOwk), Rho_EMO, sigma_EMO),
    c(p_DCS, p_DCSwk2, p_DCSwk3)[part] ~ multi_normal( c(a_DCS, b1_DCSwk, b2_DCSwk), Rho_DCS, sigma_DCS),
    c(p_AUC,p_AUCwk2, p_AUCwk3)[part] ~ multi_normal( c(a_AUC,b1_AUCwk, b2_AUCwk) , Rho_AUC , sigma_AUC),
    
    a_PSS ~ normal(.5, .1),
    b1_PSSwk ~ normal(0, .1),
    b2_PSSwk ~ normal(0, .1),
    
    a_EMO ~ normal(.5, .1),
    b1_EMOwk ~ normal(0, .1),
    b2_EMOwk ~ normal(0, .1),
    
    a_DCS ~ normal(0.5, .1),
    b1_DCSwk ~ normal(0, .1),
    b2_DCSwk ~ normal(0, .1),
    
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

# Partial Mediation Model
postP <- extract.samples(m_PartMedEC)

n_samples <- length(postP$a_MuGo)
n_obs <- length(ECmodData$MuGo)

muP_MuGo_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ECmodData$part[i]
  w2 <- ECmodData$Week2[i]
  w3 <- ECmodData$Week3[i]
  
  #stress predictors
  pss <- ECmodData$PSS[i]
  emo <- ECmodData$EmoWB[i]
  dcs <- ECmodData$DCS[i]
  auc <- ECmodData$AUC[i]
  
  for (s in 1:n_samples) {
    base <- postP$p_MuGo[s, pid]
    wk2 <- postP$p_MuGowk2[s, pid] * w2
    wk3 <- postP$p_MuGowk3[s, pid] * w3
    stress <- 0
    if (!is.na(pss)) stress <- stress + postP$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postP$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postP$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postP$b_AUC[s] * auc
    
    muP_MuGo_manual[s, i] <- base + wk2 + wk3 + stress
  }
}

muP_mean <- apply(muP_MuGo_manual, 2, mean)
muP_PI <- apply(muP_MuGo_manual, 2, PI)


# Full Mediation Model
postF <- extract.samples(m_FullMedEC)

muF_MuGo_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ECmodData$part[i]
  pss <- ECmodData$PSS[i]
  emo <- ECmodData$EmoWB[i]
  dcs <- ECmodData$DCS[i]
  auc <- ECmodData$AUC[i]
  
  for (s in 1:n_samples) {
    base <- postF$p_MuGo[s, pid]
    stress <- 0
    if (!is.na(pss)) stress <- stress + postF$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postF$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postF$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postF$b_AUC[s] * auc
    
    muF_MuGo_manual[s, i] <- base + stress
  }
}

muF_mean <- apply(muF_MuGo_manual, 2, mean)
muF_PI <- apply(muF_MuGo_manual, 2, PI)


# No Mediation Model
postN <- extract.samples(m_TimeEC)

muN_MuGo_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ECmodData$part[i]
  w2 <- ECmodData$Week2[i]
  w3 <- ECmodData$Week3[i]
  
  for (s in 1:n_samples) {
    base <- postN$p_MuGo[s, pid]
    wk2 <- postN$p_MuGowk2[s, pid] * w2
    wk3 <- postN$p_MuGowk3[s, pid] * w3
    
    muN_MuGo_manual[s, i] <- base + wk2 + wk3
  }
}

muN_mean <- apply(muN_MuGo_manual, 2, mean)
muN_PI <- apply(muN_MuGo_manual, 2, PI)


#COMBINED PIC
# Open PNG device
png("Images/Expl_MuGo_Combined_PostPredPlot.png", width = 1800, height = 600, res = 150)

# Layout and styling
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2), cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.6)

# --- Partial Mediation ---
plot(muP_mean ~ ECmodData$MuGo, col = rangi2, ylim = range(muP_PI),
     main = "Partial Mediation Model", xlab = "Observed MuGo", ylab = "Predicted MuGo")
abline(a = 0, b = 1, lty = 2)
for (i in 1:n_obs) {
  lines(rep(ECmodData$MuGo[i], 2), muP_PI[, i], col = rangi2)
}

# --- Full Mediation ---
plot(muF_mean ~ ECmodData$MuGo, col = rangi2, ylim = range(muF_PI),
     main = "Full Mediation Model", xlab = "Observed MuGo", ylab = "Predicted MuGo")
abline(a = 0, b = 1, lty = 2)
for (i in 1:n_obs) {
  lines(rep(ECmodData$MuGo[i], 2), muF_PI[, i], col = rangi2)
}

# --- No Mediation ---
plot(muN_mean ~ ECmodData$MuGo, col = rangi2, ylim = range(muN_PI),
     main = "No Mediation Model", xlab = "Observed MuGo", ylab = "Predicted MuGo")
abline(a = 0, b = 1, lty = 2)
for (i in 1:n_obs) {
  lines(rep(ECmodData$MuGo[i], 2), muN_PI[, i], col = rangi2)
}

# Close device
dev.off()


# Compare WAICs
MuGoCompare <- compare(m_PartMedEC, m_FullMedEC, m_TimeEC, func = WAIC)
MuGoComp_df <- as.data.frame(MuGoCompare)
MuGoComp_df$model <- rownames(MuGoComp_df)

# Rename models for plotting
MuGoComp_df$model[grepl("m_TimeEC", MuGoComp_df$model)] <- "No Mediation"
MuGoComp_df$model[grepl("m_PartMedEC", MuGoComp_df$model)] <- "Partial Mediation"
MuGoComp_df$model[grepl("m_FullMedEC", MuGoComp_df$model)] <- "Full Mediation"

# Move model column to front
MuGoComp_df <- MuGoComp_df[, c("model", setdiff(names(MuGoComp_df), "model"))]

# Save comparison table
write.csv(MuGoComp_df, "Output/Expl_MuGo_ModelComp.csv", row.names = FALSE)

