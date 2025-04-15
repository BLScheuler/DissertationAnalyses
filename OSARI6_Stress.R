#OSARI Step 6
#Stress mediation 

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)
library(dplyr)

#set working directory

#load exgSS modeling data
load("BeestsSubjMeans.RData")

#rename week to match, and only keep columns I'm interested in
SSRTdata <- beests_subj_means %>%
  rename(Week = week) %>%  # match the capitalization
  mutate(Week = as.integer(Week)) %>% #temporarily make it integer for left_join
  select(subject, Week, SSRT)

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
  left_join(SSRTdata, by = c("subject", "Week"))

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
  DCS = -1*(standardize(AllDat$DiurnalSlope)), #flipped to match direction of stress in model
  AUC = standardize(AllDat$AUCg), #standard normalization
  SSRT = AllDat$SSRT
)

#Set seed for reproducibility
set.seed(123)  


### ULAM MODELS

# PARTRIAL MEDIAION MODEL (TIME -> EXECUTIVE CONTROL <- STRESS; TIME -> STRESS)
m_PartMedEC <- ulam(
  alist(
    ## TIME -> EXECUTIVE CONTROL <- STRESS
    #distribution for EC parameter 
    SSRT ~ dnorm(mu_SSRT, sigma_SSRT),
    #Set up participant change in SSRT over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_SSRT <- p_SSRT[part] + p_SSRTwk2[part]*Week2 + p_SSRTwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_SSRT,p_SSRTwk2, p_SSRTwk3)[part] ~ multi_normal( c(a_SSRT,b1_SSRTwk, b2_SSRTwk) , Rho_SSRT , sigma_indSSRT),
    
    #Matzke et al 2021 for SSRT priors
    a_SSRT ~ normal(.169, .06),
    b1_SSRTwk ~ normal(0, .1),
    b2_SSRTwk ~ normal(0, .1),
    
    #Priors for stress
    b_PSS ~ normal(0, .1),
    b_EmoWB ~ normal(0, .1), 
    b_DCS ~ normal(0, .1), 
    b_AUC ~ normal(0, .1),
    
    sigma_SSRT ~ exponential(1),
    sigma_indSSRT ~ exponential(1),
    Rho_SSRT ~ lkj_corr(2),
    
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
    
  ), data = ECmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_PartMedEC, depth=3)




#FULL MEDIATION MODEL (TIME -> STRESS -> EXECUTIVE CONTROL)
m_FullMedEC <- ulam(
  alist(
    ## STRESS -> EXECUTIVE CONTROL
    #distribution for EC parameter 
    SSRT ~ dnorm(mu_SSRT, sigma_SSRT),
    #Set up participant change over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3
    mu_SSRT <- p_SSRT[part] + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    
    #Matzke et al 2021 for SSRT priors
    p_SSRT[part] ~ normal(a_SSRT, sigma_indSSRT),
    a_SSRT ~ normal(.169,.06),
    
    #Priors for stress
    b_PSS ~ normal(0,.1),
    b_EmoWB ~ normal(0,.1), 
    b_DCS ~ normal(0, .1), 
    b_AUC ~ normal(0, .1),
    
    sigma_SSRT ~ exponential(1),
    sigma_indSSRT ~ exponential(1),
    
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
    
    a_PSS ~ normal(0,.1),
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
    
  ), data = ECmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_FullMedEC, depth=3)




#NO MEDIATION MODEL (TIME -> EXECUTIVE CONTROL; TIME -> STRESS)
m_TimeEC <- ulam(
  alist(
    ## TIME -> EXECUTIVE CONTROL
    #distribution for EC parameter 
    SSRT ~ dnorm(mu_SSRT, sigma_SSRT),
    #Set up participant change over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3
    mu_SSRT <- p_SSRT[part] + p_SSRTwk2[part]*Week2 + p_SSRTwk3[part]*Week3,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_SSRT,p_SSRTwk2, p_SSRTwk3)[part] ~ multi_normal( c(a_SSRT,b1_SSRTwk, b2_SSRTwk) , Rho_SSRT , sigma_indSSRT),
    
    #Matzke et al 2021 for SSRT priors
    a_SSRT ~ normal(.169,.06),
    b1_SSRTwk ~ normal(0,.1),
    b2_SSRTwk ~ normal(0,.1),
    
    
    sigma_SSRT ~ exponential(1),
    sigma_indSSRT ~ exponential(1),
    Rho_SSRT ~ lkj_corr(2),
    
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
    
  ), data = ECmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_TimeEC, depth=3)


### POSTERIOR PREDICTIVE PLOTS

#How well did the Partial Mediation Model approximate the posterior distribution?

postP <- extract.samples(m_PartMedEC)

# Number of samples and observations
n_samples <- length(postP$a_SSRT)
n_obs <- length(ECmodData$SSRT)

# Storage for predicted means
muP_SSRT_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

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
    base <- postP$p_SSRT[s, pid]
    wk2 <- postP$p_SSRTwk2[s, pid] * w2
    wk3 <- postP$p_SSRTwk3[s, pid] * w3
    
    stress <- 0
    if (!is.na(pss)) stress <- stress + postP$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postP$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postP$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postP$b_AUC[s] * auc
    
    muP_SSRT_manual[s, i] <- base + wk2 + wk3 + stress
  }
}

muP_mean <- apply(muP_SSRT_manual, 2, mean)
muP_PI <- apply(muP_SSRT_manual, 2, PI)

#Plot it out
png("Images/SSRT_Hyp3_PartMed_PostPredPlot.png", width = 600, height = 600)
plot(muP_mean ~ECmodData$SSRT, col = rangi2, ylim = range(muP_PI),
     main= "Posterior Predictive Plots for SSRTs from Partial Mediation Model", 
     xlab = "Observed SSRT",  ylab = "Predicted SSRT")
abline(a=0, b=1, lty=2)
for(i in 1:length(ECmodData$SSRT)) {
  lines(rep(ECmodData$SSRT[i], 2), muP_PI[, i], col = rangi2)}
dev.off()

#How well did the Full Mediation Model approximate the posterior distribution?
postF <- extract.samples(m_FullMedEC)

muF_SSRT_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ECmodData$part[i]
  pss <- ECmodData$PSS[i]
  emo <- ECmodData$EmoWB[i]
  dcs <- ECmodData$DCS[i]
  auc <- ECmodData$AUC[i]
  
  for (s in 1:n_samples) {
    base <- postF$p_SSRT[s, pid]
    stress <- 0
    if (!is.na(pss)) stress <- stress + postF$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postF$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postF$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postF$b_AUC[s] * auc
    
    muF_SSRT_manual[s, i] <- base + stress
  }
}

muF_mean <- apply(muF_SSRT_manual, 2, mean)
muF_PI <- apply(muF_SSRT_manual, 2, PI)

#Plot it out
png("Images/SSRT_Hyp3_FullMed_PostPredPlot.png", width = 600, height = 600)
plot(muF_mean ~ECmodData$SSRT, col = rangi2, ylim = range(muF_PI),
     main= "Posterior Predictive Plots for SSRTs from Full Mediation Model", 
     xlab = "Observed SSRT",  ylab = "Predicted SSRT")
abline(a=0, b=1, lty=2)
for(i in 1:length(ECmodData$SSRT)) {
  lines(rep(ECmodData$SSRT[i], 2), muF_PI[, i], col = rangi2)}
dev.off()

#How well did the No Mediation Model approximate the posterior distribution?
postN <- extract.samples(m_TimeEC)

muN_SSRT_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ECmodData$part[i]
  w2 <- ECmodData$Week2[i]
  w3 <- ECmodData$Week3[i]
  
  for (s in 1:n_samples) {
    base <- postN$p_SSRT[s, pid]
    wk2 <- postN$p_SSRTwk2[s, pid] * w2
    wk3 <- postN$p_SSRTwk3[s, pid] * w3
    
    muN_SSRT_manual[s, i] <- base + wk2 + wk3
  }
}

muN_mean <- apply(muN_SSRT_manual, 2, mean)
muN_PI <- apply(muN_SSRT_manual, 2, PI)

#Plot it out
png("Images/SSRT_Hyp3_NoMed_PostPredPlot.png", width = 600, height = 600)
plot(muN_mean ~ECmodData$SSRT, col = rangi2, ylim = range(muN_PI),
     main= "Posterior Predictive Plots for SSRTs from No Mediation Model", 
     xlab = "Observed SSRT",  ylab = "Predicted SSRT")
abline(a=0, b=1, lty=2)
for(i in 1:length(ECmodData$SSRT)) {
  lines(rep(ECmodData$SSRT[i], 2), muN_PI[, i], col = rangi2)}
dev.off()


# All the plots in one image
# Open PNG device with higher resolution
png("Images/SSRT_Hyp3_Combined_PostPredPlot.png", width = 1800, height = 600, res = 150)

# Set up 1 row, 3 columns layout and increase overall text size
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2), cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.6)

# --- Partial Mediation Model ---
plot(muP_mean ~ ECmodData$SSRT, col = rangi2, ylim = range(muP_PI),
     main = "Partial Mediation Model", xlab = "Observed SSRT", ylab = "Predicted SSRT")
abline(a = 0, b = 1, lty = 2)
for (i in 1:length(ECmodData$SSRT)) {
  lines(rep(ECmodData$SSRT[i], 2), muP_PI[, i], col = rangi2)
}

# --- Full Mediation Model ---
plot(muF_mean ~ ECmodData$SSRT, col = rangi2, ylim = range(muF_PI),
     main = "Full Mediation Model", xlab = "Observed SSRT", ylab = "Predicted SSRT")
abline(a = 0, b = 1, lty = 2)
for (i in 1:length(ECmodData$SSRT)) {
  lines(rep(ECmodData$SSRT[i], 2), muF_PI[, i], col = rangi2)
}

# --- No Mediation Model ---
plot(muN_mean ~ ECmodData$SSRT, col = rangi2, ylim = range(muN_PI),
     main = "No Mediation Model", xlab = "Observed SSRT", ylab = "Predicted SSRT")
abline(a = 0, b = 1, lty = 2)
for (i in 1:length(ECmodData$SSRT)) {
  lines(rep(ECmodData$SSRT[i], 2), muN_PI[, i], col = rangi2)
}

# Close the device
dev.off()


### MODEL COMPARISON

#Compare WAICs
SSRTcompare <- compare(m_PartMedEC, m_FullMedEC, m_TimeEC, func = WAIC)
SSRTcomp_df <- as.data.frame(SSRTcompare)
SSRTcomp_df$model <- rownames(SSRTcomp_df)

#Renaming to look cleaner on plots
SSRTcomp_df$model[grepl("m_TimeEC", SSRTcomp_df$model)] <- "No Mediation"
SSRTcomp_df$model[grepl("m_PartMedEC", SSRTcomp_df$model)] <- "Partial Mediation"
SSRTcomp_df$model[grepl("m_FullMedEC", SSRTcomp_df$model)] <- "Full Mediation"

#move model column to front
SSRTcomp_df <- SSRTcomp_df[, c("model", setdiff(names(SSRTcomp_df), "model"))]

# Save to CSV
write.csv(SSRTcomp_df, "Output/Hyp3_SSRT_ModelComp.csv", row.names = FALSE)

#Rethinking Plot 
plot(compare(m_PartMedEC, m_FullMedEC, m_TimeEC, func = WAIC))  
#See Statistical Rethinking, Chapter 7, code line 7.29

#Weight Plot
SSRT_timePlot<- ggplot(SSRTcomp_df, aes(x = reorder(model, -weight), y = weight)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = round(weight, 2)), vjust = -0.5) +
  labs(title = "SSRT Model Comparison via WAIC",
       x = "Model", y = "Model Weight") +
  theme_light(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
#Save the image
ggsave("Images/Hyp3_SSRTweightPlot.png", SSRT_timePlot, width = 6, height = 6, dpi = 300)








