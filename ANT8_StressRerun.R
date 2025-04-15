#ANT Step 8
#Re run stress without poorly fit datapoint


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

#Filter out subject 5, week 1
AllDat <- AllDat[!(AllDat$Subject == "5" & ANTdata$Week == 1), ]


#Set up data for modeling
ITmodData <- list(
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

# PARTRIAL MEDIAION MODEL (TIME -> ATTENTIONAL CONTROL <- STRESS; TIME -> STRESS)
m_PartMedIT <- ulam(
  alist(
    ## TIME -> ATTENTIONAL CONTROL <- STRESS
    #distribution for IT parameter 
    IT ~ dnorm(mu_IT, sigma_IT),
    #Set up participant change in IT over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3 AND stress
    mu_IT <- p_IT[part] + p_ITwk2[part]*Week2 + p_ITwk3[part]*Week3 + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_IT, p_ITwk2, p_ITwk3)[part] ~ multi_normal( c(a_IT,b1_ITwk, b2_ITwk) , Rho_IT, sigma_indIT),
    
    #IT Priors
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
    
  ), data = ITmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_PartMedIT, depth=3)


#FULL MEDIATION MODEL (TIME -> STRESS -> ATTENTIONAL CONTROL)
m_FullMedIT <- ulam(
  alist(
    ## STRESS -> ATTENTIONAL CONTROL
    #distribution for IT parameter 
    IT ~ dnorm(mu_IT, sigma_IT),
    #Set up participant change from stress
    mu_IT <- p_IT[part] + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    #IT priors
    p_IT[part] ~ normal(a_IT, sigma_indIT),
    a_IT ~ normal(0, 1),  #Based on standardized values
    
    #Priors for stress
    b_PSS ~ normal(0, 1),
    b_EmoWB ~ normal(0, 1), 
    b_DCS ~ normal(0, 1), 
    b_AUC ~ normal(0, 1), 
    
    sigma_IT ~ exponential(1),
    sigma_indIT ~ exponential(1),
    
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
    
  ), data = ITmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_FullMedIT, depth=3)



#NO MEDIATION MODEL (TIME -> ATTENTIONAL CONTROL; TIME -> STRESS)
m_TimeIT <- ulam(
  alist(
    ## TIME -> ATTENTIONAL CONTROL
    #distribution for IT parameter 
    IT ~ dnorm(mu_IT, sigma_IT),
    #Set up participant change over time
    #Week 1 (where wk2 and wk3 are zero) plus changes for wk 2 and wk 3
    mu_IT <- p_IT[part] + p_ITwk2[part]*Week2 + p_ITwk3[part]*Week3,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_IT, p_ITwk2, p_ITwk3)[part] ~ multi_normal( c(a_IT, b1_ITwk, b2_ITwk) , Rho_IT, sigma_indIT),
    
    #IT priors
    a_IT ~ normal(0, 1),
    b1_ITwk ~ normal(0, 1),
    b2_ITwk ~ normal(0, 1),
    
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
    
    c(p_PSS, p_PSSwk2, p_PSSwk3)[part] ~ multi_normal( c(a_PSS, b1_PSSwk, b2_PSSwk) , Rho_PSS , sigma_PSS),
    c(p_EmoWB, p_EMOwk2, p_EMOwk3)[part] ~ multi_normal( c(a_EMO, b1_EMOwk, b2_EMOwk) , Rho_EMO , sigma_EMO),
    c(p_DCS, p_DCSwk2, p_DCSwk3)[part] ~ multi_normal( c(a_DCS, b1_DCSwk, b2_DCSwk) , Rho_DCS , sigma_DCS), 
    c(p_AUC, p_AUCwk2, p_AUCwk3)[part] ~ multi_normal( c(a_AUC, b1_AUCwk, b2_AUCwk) , Rho_AUC , sigma_AUC), 
    
    a_PSS ~ normal(.5, .1),
    b1_PSSwk ~ normal(0, .1),
    b2_PSSwk ~ normal(0, .1),
    
    a_EMO ~ normal(.5, .1),
    b1_EMOwk ~ normal(0, .1),
    b2_EMOwk ~ normal(0, .1),
    
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
    
  ), data = ITmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_TimeIT, depth=3)

### POSTERIOR PREDICTIVE PLOTS

#How well did the Partial Mediation Model approximate the posterior distribution?


#Extract posterior samples for Partial Mediation model
postP <- extract.samples(m_PartMedIT)

# Number of samples and observations
n_samples <- length(postP$a_IT)
n_obs <- length(ITmodData$IT)

# Storage for predicted means
muP_IT_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ITmodData$part[i]
  w2 <- ITmodData$Week2[i]
  w3 <- ITmodData$Week3[i]
  
  pss <- ITmodData$PSS[i]
  emo <- ITmodData$EmoWB[i]
  dcs <- ITmodData$DCS[i]
  auc <- ITmodData$AUC[i]
  
  for (s in 1:n_samples) {
    base <- postP$p_IT[s, pid]
    wk2 <- postP$p_ITwk2[s, pid] * w2
    wk3 <- postP$p_ITwk3[s, pid] * w3
    
    stress <- 0
    if (!is.na(pss)) stress <- stress + postP$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postP$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postP$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postP$b_AUC[s] * auc
    
    muP_IT_manual[s, i] <- base + wk2 + wk3 + stress
  }
}

muP_mean <- apply(muP_IT_manual, 2, mean)
muP_PI <- apply(muP_IT_manual, 2, PI)

# #Plot
# png("Images/IT_Hyp3_PartMed_PostPredPlot.png", width = 600, height = 600)
# plot(muP_mean ~ ITmodData$IT, col = rangi2, ylim = range(muP_PI),
#      main = "Posterior Predictive Plots for IT from Partial Mediation Model", 
#      xlab = "Observed IT", ylab = "Predicted IT")
# abline(a = 0, b = 1, lty = 2)
# for (i in 1:length(ITmodData$IT)) {
#   lines(rep(ITmodData$IT[i], 2), muP_PI[, i], col = rangi2)
# }
# dev.off()


#How well did the Full Mediation Model approximate the posterior distribution?
postF <- extract.samples(m_FullMedIT)
muF_IT_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ITmodData$part[i]
  pss <- ITmodData$PSS[i]
  emo <- ITmodData$EmoWB[i]
  dcs <- ITmodData$DCS[i]
  auc <- ITmodData$AUC[i]
  
  for (s in 1:n_samples) {
    base <- postF$p_IT[s, pid]
    stress <- 0
    if (!is.na(pss)) stress <- stress + postF$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postF$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postF$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postF$b_AUC[s] * auc
    
    muF_IT_manual[s, i] <- base + stress
  }
}

muF_mean <- apply(muF_IT_manual, 2, mean)
muF_PI <- apply(muF_IT_manual, 2, PI)

# png("Images/IT_Hyp3_FullMed_PostPredPlot.png", width = 600, height = 600)
# plot(muF_mean ~ ITmodData$IT, col = rangi2, ylim = range(muF_PI),
#      main = "Posterior Predictive Plots for IT from Full Mediation Model", 
#      xlab = "Observed IT", ylab = "Predicted IT")
# abline(a = 0, b = 1, lty = 2)
# for (i in 1:length(ITmodData$IT)) {
#   lines(rep(ITmodData$IT[i], 2), muF_PI[, i], col = rangi2)
# }
# dev.off()



#How well did the No Mediation Model approximate the posterior distribution?
postN <- extract.samples(m_TimeIT)
muN_IT_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- ITmodData$part[i]
  w2 <- ITmodData$Week2[i]
  w3 <- ITmodData$Week3[i]
  
  for (s in 1:n_samples) {
    base <- postN$p_IT[s, pid]
    wk2 <- postN$p_ITwk2[s, pid] * w2
    wk3 <- postN$p_ITwk3[s, pid] * w3
    muN_IT_manual[s, i] <- base + wk2 + wk3
  }
}

muN_mean <- apply(muN_IT_manual, 2, mean)
muN_PI <- apply(muN_IT_manual, 2, PI)



# ## I want all the plots in one image
# png("Images/IT_Hyp3_NoMed_PostPredPlot.png", width = 600, height = 600)
# plot(muN_mean ~ ITmodData$IT, col = rangi2, ylim = range(muN_PI),
#      main = "Posterior Predictive Plots for IT from No Mediation Model", 
#      xlab = "Observed IT", ylab = "Predicted IT")
# abline(a = 0, b = 1, lty = 2)
# for (i in 1:length(ITmodData$IT)) {
#   lines(rep(ITmodData$IT[i], 2), muN_PI[, i], col = rangi2)
# }
# dev.off()
# 
# png("Images/IT_Hyp3_Combined_PostPredPlot.png", width = 1800, height = 600, res = 150)
# par(mfrow = c(1, 3), mar = c(5, 5, 4, 2), cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.6)
# 
# # Partial
# plot(muP_mean ~ ITmodData$IT, col = rangi2, ylim = range(muP_PI),
#      main = "Partial Mediation Model", xlab = "Observed IT", ylab = "Predicted IT")
# abline(a = 0, b = 1, lty = 2)
# for (i in 1:length(ITmodData$IT)) {
#   lines(rep(ITmodData$IT[i], 2), muP_PI[, i], col = rangi2)
# }
# 
# # Full
# plot(muF_mean ~ ITmodData$IT, col = rangi2, ylim = range(muF_PI),
#      main = "Full Mediation Model", xlab = "Observed IT", ylab = "Predicted IT")
# abline(a = 0, b = 1, lty = 2)
# for (i in 1:length(ITmodData$IT)) {
#   lines(rep(ITmodData$IT[i], 2), muF_PI[, i], col = rangi2)
# }
# 
# # No
# plot(muN_mean ~ ITmodData$IT, col = rangi2, ylim = range(muN_PI),
#      main = "No Mediation Model", xlab = "Observed IT", ylab = "Predicted IT")
# abline(a = 0, b = 1, lty = 2)
# for (i in 1:length(ITmodData$IT)) {
#   lines(rep(ITmodData$IT[i], 2), muN_PI[, i], col = rangi2)
# }
# dev.off()



### MODEL COMPARISON

#Compare WAICs
ITcompare <- compare(m_PartMedIT, m_FullMedIT, m_TimeIT, func = WAIC)
ITcomp_df <- as.data.frame(ITcompare)
ITcomp_df$model <- rownames(ITcomp_df)

#Renaming to look cleaner on plots
ITcomp_df$model[grepl("m_TimeIT", ITcomp_df$model)] <- "No Mediation"
ITcomp_df$model[grepl("m_PartMedIT", ITcomp_df$model)] <- "Partial Mediation"
ITcomp_df$model[grepl("m_FullMedIT", ITcomp_df$model)] <- "Full Mediation"

#Move model column to front
ITcomp_df <- ITcomp_df[, c("model", setdiff(names(ITcomp_df), "model"))]

# Save to CSV
write.csv(ITcomp_df, "Output/Hyp3Check_IT_ModelComp.csv", row.names = FALSE)

#Rethinking Plot 
#plot(compare(m_PartMedIT, m_FullMedIT, m_TimeIT, func = WAIC))  
#See Statistical Rethinking, Chapter 7, code line 7.29

# #Weight Plot
# IT_weightPlot <- ggplot(ITcomp_df, aes(x = reorder(model, -weight), y = weight)) +
#   geom_bar(stat = "identity", fill = "skyblue") +
#   geom_text(aes(label = round(weight, 2)), vjust = -0.5) +
#   labs(title = "IT Model Comparison via WAIC",
#        x = "Model", y = "Model Weight") +
#   theme_light(base_size = 14) +
#   theme(
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank()
#   )
# 
# #Save the image
# ggsave("Images/Hyp3_ITweightPlot.png", IT_weightPlot, width = 6, height = 6, dpi = 300)
