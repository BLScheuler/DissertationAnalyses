#DNB Step 6
#Stress mediation


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
WCmodData <- list(
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

# PARTRIAL MEDIAION MODEL (TIME -> WORKLOAD CAPACITY <- STRESS; TIME -> STRESS)
m_PartMedWC <- ulam(
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
    
  ), data = WCmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_PartMedWC, depth=3)


#FULL MEDIATION MODEL (TIME -> STRESS -> WORKLOAD CAPACITY)
m_FullMedWC <- ulam(
  alist(
    ## STRESS -> WORKLOAD CAPACITY
    #distribution for WC parameter 
    WC ~ dnorm(mu_WC, sigma_WC),
    #Set up participant change from stress
    mu_WC <- p_WC[part] + b_PSS*PSS + b_EmoWB*EmoWB + b_DCS*DCS + b_AUC*AUC,
    
    
    #WC priors
    p_WC[part] ~ normal(a_WC, sigma_indWC),
    a_WC ~ normal(3, 0.75),
    
    #Priors for stress
    b_PSS ~ normal(0,.1),
    b_EmoWB ~ normal(0,.1), 
    b_DCS ~ normal(0, .1), 
    b_AUC ~ normal(0, .1), 
    
    sigma_WC ~ exponential(1),
    sigma_indWC ~ exponential(1),
    
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
    
  ), data = WCmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_FullMedWC, depth=3)




#NO MEDIATION MODEL (TIME -> WORKLOAD CAPACITY; TIME -> STRESS)
m_TimeWC <- ulam(
  alist(
    ## TIME -> WORKLOAD CAPACITY
    #distribution for WC parameter 
    WC ~ dnorm(mu_WC, sigma_WC),
    #Set up participant change over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3
    mu_WC <- p_WC[part] + p_WCwk2[part]*Week2 + p_WCwk3[part]*Week3,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_WC,p_WCwk2, p_WCwk3)[part] ~ multi_normal( c(a_WC,b1_WCwk, b2_WCwk) , Rho_WC, sigma_indWC),
    
    #WC priors
    a_WC ~ normal(3, 0.75),
    b1_WCwk ~ normal(0, 0.5),
    b2_WCwk ~ normal(0, 0.5),
    
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
    
  ), data = WCmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_TimeWC, depth=3)




### POSTERIOR PREDICTIVE PLOTS

#How well did the Partial Mediation Model approximate the posterior distribution?

postP <- extract.samples(m_PartMedWC)

# Number of samples and observations
n_samples <- length(postP$a_WC)
n_obs <- length(WCmodData$WC)

# Storage for predicted means
muP_WC_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- WCmodData$part[i]
  w2 <- WCmodData$Week2[i]
  w3 <- WCmodData$Week3[i]
  
  # Pull stress predictors (may be NA)
  pss <- WCmodData$PSS[i]
  emo <- WCmodData$EmoWB[i]
  dcs <- WCmodData$DCS[i]
  auc <- WCmodData$AUC[i]
  
  for (s in 1:n_samples) {
    base <- postP$p_WC[s, pid]
    wk2 <- postP$p_WCwk2[s, pid] * w2
    wk3 <- postP$p_WCwk3[s, pid] * w3
    
    stress <- 0
    if (!is.na(pss)) stress <- stress + postP$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postP$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postP$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postP$b_AUC[s] * auc
    
    muP_WC_manual[s, i] <- base + wk2 + wk3 + stress
  }
}

muP_mean <- apply(muP_WC_manual, 2, mean)
muP_PI <- apply(muP_WC_manual, 2, PI)

#Plot it out
png("Images/WC_Hyp3_PartMed_PostPredPlot.png", width = 600, height = 600)
plot(muP_mean ~WCmodData$WC, col = rangi2, ylim = range(muP_PI),
     main= "Posterior Predictive Plots for Workload Capacity from Partial Mediation Model", 
     xlab = "Observed WC",  ylab = "Predicted WC")
abline(a=0, b=1, lty=2)
for(i in 1:length(WCmodData$WC)) {
  lines(rep(WCmodData$WC[i], 2), muP_PI[, i], col = rangi2)}
dev.off()

#How well did the Full Mediation Model approximate the posterior distribution?
postF <- extract.samples(m_FullMedWC)

muF_WC_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- WCmodData$part[i]
  pss <- WCmodData$PSS[i]
  emo <- WCmodData$EmoWB[i]
  dcs <- WCmodData$DCS[i]
  auc <- WCmodData$AUC[i]
  
  for (s in 1:n_samples) {
    base <- postF$p_WC[s, pid]
    stress <- 0
    if (!is.na(pss)) stress <- stress + postF$b_PSS[s] * pss
    if (!is.na(emo)) stress <- stress + postF$b_EmoWB[s] * emo
    if (!is.na(dcs)) stress <- stress + postF$b_DCS[s] * dcs
    if (!is.na(auc)) stress <- stress + postF$b_AUC[s] * auc
    
    muF_WC_manual[s, i] <- base + stress
  }
}

muF_mean <- apply(muF_WC_manual, 2, mean)
muF_PI <- apply(muF_WC_manual, 2, PI)

#Plot it out
png("Images/WC_Hyp3_FullMed_PostPredPlot.png", width = 600, height = 600)
plot(muF_mean ~WCmodData$WC, col = rangi2, ylim = range(muF_PI),
     main= "Posterior Predictive Plots for Workload Capacity from Full Mediation Model", 
     xlab = "Observed WC",  ylab = "Predicted WC")
abline(a=0, b=1, lty=2)
for(i in 1:length(WCmodData$WC)) {
  lines(rep(WCmodData$WC[i], 2), muF_PI[, i], col = rangi2)}
dev.off()

#How well did the No Mediation Model approximate the posterior distribution?
postN <- extract.samples(m_TimeWC)

muN_WC_manual <- matrix(NA, nrow = n_samples, ncol = n_obs)

for (i in 1:n_obs) {
  pid <- WCmodData$part[i]
  w2 <- WCmodData$Week2[i]
  w3 <- WCmodData$Week3[i]
  
  for (s in 1:n_samples) {
    base <- postN$p_WC[s, pid]
    wk2 <- postN$p_WCwk2[s, pid] * w2
    wk3 <- postN$p_WCwk3[s, pid] * w3
    
    muN_WC_manual[s, i] <- base + wk2 + wk3
  }
}

muN_mean <- apply(muN_WC_manual, 2, mean)
muN_PI <- apply(muN_WC_manual, 2, PI)


#Plot it out
png("Images/WC_Hyp3_NoMed_PostPredPlot.png", width = 600, height = 600)
plot(muN_mean ~WCmodData$WC, col = rangi2, ylim = range(muN_PI),
     main= "Posterior Predictive Plots for Workload Capacity from No Mediation Model", 
     xlab = "Observed WC",  ylab = "Predicted WC")
abline(a=0, b=1, lty=2)
for(i in 1:length(WCmodData$WC)) {
  lines(rep(WCmodData$WC[i], 2), muN_PI[, i], col = rangi2)}
dev.off()


##all the plots in one image
# Open PNG device with higher resolution
png("Images/WC_Hyp3_Combined_PostPredPlot.png", width = 1800, height = 600, res = 150)

# Set up 1 row, 3 columns layout and increase overall text size
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2), cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.6)

# --- Partial Mediation Model ---
plot(muP_mean ~ WCmodData$WC, col = rangi2, ylim = range(muP_PI),
     main = "Partial Mediation Model", xlab = "Observed WC", ylab = "Predicted WC")
abline(a = 0, b = 1, lty = 2)
for (i in 1:length(WCmodData$WC)) {
  lines(rep(WCmodData$WC[i], 2), muP_PI[, i], col = rangi2)
}

# --- Full Mediation Model ---
plot(muF_mean ~ WCmodData$WC, col = rangi2, ylim = range(muF_PI),
     main = "Full Mediation Model", xlab = "Observed WC", ylab = "Predicted WC")
abline(a = 0, b = 1, lty = 2)
for (i in 1:length(WCmodData$WC)) {
  lines(rep(WCmodData$WC[i], 2), muF_PI[, i], col = rangi2)
}

# --- No Mediation Model ---
plot(muN_mean ~ WCmodData$WC, col = rangi2, ylim = range(muN_PI),
     main = "No Mediation Model", xlab = "Observed WC", ylab = "Predicted WC")
abline(a = 0, b = 1, lty = 2)
for (i in 1:length(WCmodData$WC)) {
  lines(rep(WCmodData$WC[i], 2), muN_PI[, i], col = rangi2)
}

# Close the device
dev.off()


### MODEL COMPARISON

#Compare WAICs
WCcompare <- compare(m_PartMedWC, m_FullMedWC, m_TimeWC, func = WAIC)
WCcomp_df <- as.data.frame(WCcompare)
WCcomp_df$model <- rownames(WCcomp_df)

#Renaming to look cleaner on plots
WCcomp_df$model[grepl("m_TimeWC", WCcomp_df$model)] <- "No Mediation"
WCcomp_df$model[grepl("m_PartMedWC", WCcomp_df$model)] <- "Partial Mediation"
WCcomp_df$model[grepl("m_FullMedWC", WCcomp_df$model)] <- "Full Mediation"

#move model column to front
WCcomp_df <- WCcomp_df[, c("model", setdiff(names(WCcomp_df), "model"))]

# Save to CSV
write.csv(WCcomp_df, "Output/Hyp3_WC_ModelComp.csv", row.names = FALSE)

#Rethinking Plot 
plot(compare(m_PartMedWC, m_FullMedWC, m_TimeWC, func = WAIC))  
#See Statistical Rethinking, Chapter 7, code line 7.29

#Weight Plot
WC_weightPlot<- ggplot(WCcomp_df, aes(x = reorder(model, -weight), y = weight)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = round(weight, 2)), vjust = -0.5) +
  labs(title = "WC Model Comparison via WAIC",
       x = "Model", y = "Model Weight") +
  theme_light(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
#Save the image
ggsave("Images/Hyp3_WCweightPlot.png", WC_weightPlot, width = 6, height = 6, dpi = 300)

