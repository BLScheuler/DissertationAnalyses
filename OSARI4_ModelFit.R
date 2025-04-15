#OSARI Step 4
#Test model fit

#set working directory

library(coda)     # For computing HPD intervals (HDI)
library(dplyr)    # For data manipulation
library(ggplot2)  # For the violin plots
library(dplyr)    # For data manipulation
library(tidyr)    # For the correlations
library(purrr)    # For the correlations


# Load DMC material
source("dmc/dmc.R")
source("dmc/dmc_model.R")
source("dmc/dmc_sampling.R")
source("dmc/dmc_hierarchical.R")
source("dmc/dmc_analysis.R")

# Load the EXG-SS (BEESTS) model
load_model("EXG-SS", "exgSS.R")


# Load cleaned OSARI data
OSARIdata <- read.csv("CleanOSARI.csv")

#Load posterior parameters from exgSS modeling
load("osari_posterior.RData")   # This loads samples_subject_week
OSARIpost <- samples_subject_week

#turn into dataframe
post_df <- c()
for (sj in 1:5){
  for (wk in 1:3) {
    mu.true <- c(OSARIpost[[sj]][[wk]]$theta[,"mu.true",])
    sigma.true <- c(OSARIpost[[sj]][[wk]]$theta[,"sigma.true",])
    tau.true <- c(OSARIpost[[sj]][[wk]]$theta[,"tau.true",])
    muS <- c(OSARIpost[[sj]][[wk]]$theta[,"muS",])
    sigmaS <- c(OSARIpost[[sj]][[wk]]$theta[,"sigmaS",])
    tauS <- c(OSARIpost[[sj]][[wk]]$theta[,"tauS",])
    
    s1post <- data.frame(subject=sj, week=wk, mu.true=mu.true, sigma.true=sigma.true, tau.true, muS, sigmaS, tauS)
    s1post_ix <- sample(dim(s1post)[1], 500, replace=FALSE)
    s1post <- s1post[s1post_ix,]
    post_df <- rbind(post_df, s1post)
  }}

#Turn subject and week into factors
post_df$subject <- factor(post_df$subject)
post_df$week <- factor(post_df$week, ordered=TRUE)

#Reproducibility for simulations
set.seed(123)

# Define the exgSS model (from modeling code and dmc page)
model <- model.dmc(
  # SS stands for trial type (GO or Stop-signal [SS]):
  factors=list(S= "s1", SS=c("GO","SS")),
  # NR stands for "No response", i.e., go omission & successful inhibitions:
  responses=c("NR","r1"),
  # Match scores correct responses for each GO stimulus as usual, and scores 
  # the "correct" stimulus corresponding to an NR response, but the latter has
  match.map=list(M=list(s1="r1",s1="NR")), #same s1 twice
  p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
  # No errors:
  constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001, tf=0, gf=0),
  type="exgss")

#get median for all OSARI parameters (can use them in p vector, but less appropriate)
mediParms <- post_df %>%
  group_by(subject, week) %>%
  summarise(mu = median(mu.true), sigma = median(sigma.true), tau = median(tau.true), muS = median(muS), sigmaS = median(sigmaS), tauS = median(tauS))

write.csv(mediParms, file = "Output/OSARI_MediParameters.csv", row.names = FALSE)


ppc_df <- c()
for (sj in 1:5){
  for (wk in 1:3){
    for (samp in 1:500) {
      subdat <- post_df[post_df$subject==sj & post_df$week==wk,]
      p.vector <- as.numeric(subdat[samp,c("mu.true", "sigma.true", "tau.true", "muS", "sigmaS", "tauS")])
      names(p.vector) <- c("mu.true", "sigma.true", "tau.true", "muS", "sigmaS", "tauS")
      
      check.p.vector(p.vector,model) # if works, says nothing
      data_fixed <- data.model.dmc(simulate.dmc(p.vector,model,n=1,SSD=c(Inf,.50)),model) # Check experiment code for actual stop signal time (.5)
      s1ppc_df <- data.frame(subject=sj, week=wk, trial_type=data_fixed$SS, inhibit=data_fixed$R, RT=data_fixed$RT)
      ppc_df <- rbind(ppc_df, s1ppc_df)
    }}}




### OVERALL HDIs


# Prepare HDI results dataframe for GO and failed STOP trials 
OSARI_HDIcheck_GO <- data.frame(subject = numeric(),
                                week = numeric(),
                                pct_in_HDI = numeric())

OSARI_HDIcheck_fSTOP <- data.frame(subject = numeric(),
                                   week = numeric(),
                                   pct_in_HDI = numeric())

#Need to split between GO and STOP trials, both in observed data and in simulations


#Observed Data GO Trials
OSARIobsGo <- OSARIdata %>%
  filter(trial_type == "GO")

#Observed Data Unsuccessful STOP Trials
OSARIobsStop <- OSARIdata %>%
  filter(trial_type == "STOP", inhibited == 0)

#Observed Data GO Trials
exgSSsimGo <- ppc_df %>%
  filter(trial_type == "GO")

#Observed Data Unsuccessful STOP Trials
exgSSsimStop <- ppc_df %>%
  filter(trial_type == "SS", inhibit == "r1")


#GO Trial HDI loop
for (sj in 1:5){
  for (wk in 1:3){
    #sj<- 1
    #wk <- 1
    #pull in the RTs (using dplyr)
    GOobsRTs <- OSARIobsGo %>% filter (subject==sj, week==wk) %>% pull(rt)
    GOsimRTs <- exgSSsimGo %>% filter (subject==sj, week==wk) %>% pull(RT)
    
    #Get 66% HDI of simulated RTs
    GOmcmc_sim <- as.mcmc(GOsimRTs)  #coda package needs RTs as MCMC object before it can do HPD interval
    GO_HDIbounds <- HPDinterval(GOmcmc_sim, prob = 0.66)
    GOlower <- GO_HDIbounds[1]
    GOupper <- GO_HDIbounds[2]
    
    # What percent of observed RTs fall within the HDI
    GOin_hdi <- GOobsRTs >= GOlower & GOobsRTs <= GOupper
    pct <- mean(GOin_hdi, na.rm=TRUE) * 100
    
    # Save results
    OSARI_HDIcheck_GO <- rbind(OSARI_HDIcheck_GO, c(sj, wk, pct))
  }}
names(OSARI_HDIcheck_GO) <- c("Subject", "Week", "Percent")

write.csv(OSARI_HDIcheck_GO, file = "Output/GO_HDIresults.csv", row.names = FALSE)



#STOP Trial HDI loop
for (sj in 1:5){
  for (wk in 1:3){
    
    #pull in the RTs (using dplyr)
    STOPobsRTs <- OSARIobsStop %>% filter (subject==sj, week==wk) %>% pull(rt)
    STOPsimRTs <- exgSSsimStop %>% filter (subject==sj, week==wk) %>% pull(RT)
    
    #Get 66% HDI of simulated RTs
    STOPmcmc_sim <- as.mcmc(STOPsimRTs)  #coda package needs RTs as MCMC object before it can do HPD interval
    STOP_HDIbounds <- HPDinterval(STOPmcmc_sim, prob = 0.66)
    STOPlower <- STOP_HDIbounds[1]
    STOPupper <- STOP_HDIbounds[2]
    
    # What percent of observed RTs fall within the HDI
    STOPin_hdi <- STOPobsRTs >= STOPlower & STOPobsRTs <= STOPupper
    pct <- mean(STOPin_hdi, , na.rm=TRUE) * 100
    
    # Save results
    OSARI_HDIcheck_fSTOP <- rbind(OSARI_HDIcheck_fSTOP, c(sj, wk, pct))
  }}
names(OSARI_HDIcheck_fSTOP) <- c("Subject", "Week", "Percent")
write.csv(OSARI_HDIcheck_fSTOP, file = "Output/STOP_HDIresults.csv", row.names = FALSE)


### QUANTILE HDIs

#Set quantile levels and create data frames
quantile_levels <- c(0.1, 0.3, 0.5, 0.7, 0.9)
OSARI_GO_QuantResults <- data.frame()  
OSARI_STOP_QuantResults <- data.frame()  


#GO Trial Quant Loop  
#Loops for subject and week
for (sj in 1:5){
  for (wk in 1:3){
    
    #Observed trials
    #Pull in observed GO RTs
    Qobs_GO_rts <- OSARIobsGo %>% filter (subject==sj, week==wk) %>% pull(rt)
    n_obs <- length(Qobs_GO_rts)
    
    #Simulation trials
    #Pull in simulated rts
    Qsim_GO_rts <-exgSSsimGo %>% filter(subject == sj, week == wk) %>% pull(RT)
    #Split sims into rows the same size as in our observation
    n_sim_sets <- floor(length(Qsim_GO_rts)/ n_obs) #floor keeps it from having a partial last row
    #Create a matrix for thouse sims
    QGOsim_matrix <- matrix(Qsim_GO_rts[1:n_sim_sets*n_obs], nrow = n_sim_sets, byrow = TRUE)
    
    #Quantile setup
    #Split observed trials into quantiles
    GOobs_quantiles <- quantile(Qobs_GO_rts, probs = quantile_levels, na.rm = TRUE)
    #Compute quantiles for each set of sims
    GOsim_quantiles <- t(apply(QGOsim_matrix, 1, function(x) quantile(x, probs = quantile_levels)))
    #Drop any rows with NA
    GOsim_quantiles <- GOsim_quantiles[complete.cases(GOsim_quantiles),]
    
    #HDIs
    #Compute 66% HDI for each quantile
    QGO_HDIs <- apply(GOsim_quantiles, 2, function(qvals) HPDinterval(as.mcmc(qvals), prob = 0.66))
    #Check if the obs quantiles fall in the HDIs for each quantile
    QGO_within_HDI <- mapply(function(obs, bounds) obs >= bounds[1] & obs <= bounds[2],
                             obs = GOobs_quantiles, 
                             bounds = as.data.frame(QGO_HDIs))
    
    #Store all the results in a data frame
    GOresult_row <- data.frame(
      subject = sj,
      week = wk,
      q10_obs = GOobs_quantiles[1], # obs GO rt at q10
      q10_HDI_lower = QGO_HDIs[1,1],  # lower HDI from sim for q10
      q10_HDI_upper = QGO_HDIs[2,1],  # upper HDI from sim for q10
      q10_in_HDI = QGO_within_HDI[1], # whether or not obs was in HDI
      q30_obs = GOobs_quantiles[2],
      q30_HDI_lower = QGO_HDIs[1,2],
      q30_HDI_upper = QGO_HDIs[2,2],
      q30_in_HDI = QGO_within_HDI[2],
      q50_obs = GOobs_quantiles[3],
      q50_HDI_lower = QGO_HDIs[1,3],
      q50_HDI_upper = QGO_HDIs[2,3],
      q50_in_HDI = QGO_within_HDI[3],
      q70_obs = GOobs_quantiles[4],
      q70_HDI_lower = QGO_HDIs[1,4],
      q70_HDI_upper = QGO_HDIs[2,4],
      q70_in_HDI = QGO_within_HDI[4],
      q90_obs = GOobs_quantiles[5],
      q90_HDI_lower = QGO_HDIs[1,5],
      q90_HDI_upper = QGO_HDIs[2,5],
      q90_in_HDI = QGO_within_HDI[5],
      num_in_HDI = sum(QGO_within_HDI))
    
    #Store results in data frame
    OSARI_GO_QuantResults <- rbind(OSARI_GO_QuantResults, GOresult_row)
  }
}

#Save results as csv
write.csv(OSARI_GO_QuantResults, "Output/GO_HDIquantiles.csv", row.names = FALSE)



#STOP Trial Quant Loop  
#Loops for subject and week
for (sj in 1:5){
  for (wk in 1:3){
    
    #Observed trials
    #Pull in observed STOP RTs
    Qobs_STOP_rts <- OSARIobsStop %>% filter (subject==sj, week==wk) %>% pull(rt)
    n_obs <- length(Qobs_STOP_rts)
    
    #Simulation trials
    #Pull in simulated rts
    Qsim_STOP_rts <-exgSSsimStop %>% filter(subject == sj, week == wk) %>% pull(RT)
    #Split sims into rows the same size as in our observation
    n_sim_sets <- floor(length(Qsim_STOP_rts)/ n_obs) #floor keeps it from having a partial last row
    #Create a matrix for thouse sims
    QSTOPsim_matrix <- matrix(Qsim_STOP_rts[1:n_sim_sets*n_obs], nrow = n_sim_sets, byrow = TRUE)
    
    #Quantile setup
    #Split observed trials into quantiles
    STOPobs_quantiles <- quantile(Qobs_STOP_rts, probs = quantile_levels, na.rm = TRUE)
    #Compute quantiles for each set of sims
    STOPsim_quantiles <- t(apply(QSTOPsim_matrix, 1, function(x) quantile(x, probs = quantile_levels)))
    #Drop any rows with NA
    STOPsim_quantiles <- STOPsim_quantiles[complete.cases(STOPsim_quantiles),]
    #if (nrow(GOsim_quantiles) <50) next #quality check
    
    #HDIs
    #Compute 66% HDI for each quantile
    QSTOP_HDIs <- apply(STOPsim_quantiles, 2, function(qvals) HPDinterval(as.mcmc(qvals), prob = 0.66))
    #Check if the obs quantiles fall in the HDIs for each quantile
    QSTOP_within_HDI <- mapply(function(obs, bounds) obs >= bounds[1] & obs <= bounds[2],
                               obs = STOPobs_quantiles, 
                               bounds = as.data.frame(QSTOP_HDIs))
    
    #Store all the results in a data frame
    # Store full results per quantile
    STOPresult_row <- data.frame(
      subject = sj,
      week = wk,
      q10_obs = STOPobs_quantiles[1], # obs GO rt at q10
      q10_HDI_lower = QSTOP_HDIs[1,1],  # lower HDI from sim for q10
      q10_HDI_upper = QSTOP_HDIs[2,1],  # upper HDI from sim for q10
      q10_in_HDI = QSTOP_within_HDI[1], # whether or not obs was in HDI
      q30_obs = STOPobs_quantiles[2],
      q30_HDI_lower = QSTOP_HDIs[1,2],
      q30_HDI_upper = QSTOP_HDIs[2,2],
      q30_in_HDI = QSTOP_within_HDI[2],
      q50_obs = STOPobs_quantiles[3],
      q50_HDI_lower = QSTOP_HDIs[1,3],
      q50_HDI_upper = QSTOP_HDIs[2,3],
      q50_in_HDI = QSTOP_within_HDI[3],
      q70_obs = STOPobs_quantiles[4],
      q70_HDI_lower = QSTOP_HDIs[1,4],
      q70_HDI_upper = QSTOP_HDIs[2,4],
      q70_in_HDI = QSTOP_within_HDI[4],
      q90_obs = STOPobs_quantiles[5],
      q90_HDI_lower = QSTOP_HDIs[1,5],
      q90_HDI_upper = QSTOP_HDIs[2,5],
      q90_in_HDI = QSTOP_within_HDI[5],
      num_in_HDI = sum(QSTOP_within_HDI))
    
    #Store results in data frame
    OSARI_STOP_QuantResults <- rbind(OSARI_STOP_QuantResults, STOPresult_row)
  }
}

#Save results as csv
write.csv(OSARI_STOP_QuantResults, "Output/STOP_HDIquantiles.csv", row.names = FALSE)


##PLOTS
#GO Trials
# See what the actual RTs looked like
# First, calculate subject-level mean RTs per week (for lines)
OSARIGOdata_summary <- OSARIobsGo %>%
  group_by(subject, week) %>%
  summarise(mean_rt = mean(rt), .groups = "drop")


#GO Trials - Plot Observed RTs
OGOrt_plot <- ggplot(OSARIobsGo, aes(x = factor(week), y = rt)) +
  geom_violin(fill = "lightgray", color = "black", alpha = 0.6, trim = FALSE) +
  geom_jitter(aes(color = factor(subject)), width = 0.1, size = 1.2, alpha = 0.6) +
  labs(title = "OSARI Observed GO Trial RT Distributions",
       x = "Week",
       y = "Observed Response Time (s)",
       color = "Subject") +
  coord_cartesian(ylim = c(0, NA)) +  # allows auto-scaling for upper bound
  theme_light()
# Save the file 
ggsave("Images/OGOrt_plot.png", OGOrt_plot, width = 5, height = 5, dpi = 300)

#STOP Trial- Plot Observed RTs
OSTOPrt_plot <- ggplot(OSARIobsStop, aes(x = factor(week), y = rt)) +
  geom_violin(fill = "lightgray", color = "black", alpha = 0.6, trim = FALSE) +
  geom_jitter(aes(color = factor(subject)), width = 0.1, size = 1.2, alpha = 0.6) +
  labs(title = "OSARI Observed STOP Trial RT Distributions",
       x = "Week",
       y = "Observed Response Time (s)",
       color = "Subject") +
  #coord_cartesian(ylim = c(0, NA)) +  # allows auto-scaling for upper bound
  theme_light()

# Save the file 
ggsave("Images/OSTOPrt_plot.png", OSTOPrt_plot, width = 5, height = 5, dpi = 300)



# Create clean subset for observed GO trials
obs_go_subset <- OSARIobsGo %>%
  select(subject, week, rt) %>%
  mutate(Source = "Observed")

# Create clean subset for simulated GO trials
sim_go_subset <- exgSSsimGo %>%
  select(subject, week, RT) %>%
  rename(rt = RT) %>%
  mutate(Source = "Simulated")

# Combine both into one comparison dataset
OGOcomparison_data <- rbind(obs_go_subset, sim_go_subset)

# Overlay violin plot where Simulated = violin, Observed = jitter
# Boxes split by subject, with 3 weeks inside
OGOcomparison_plot1 <- ggplot() +
  geom_violin(data = OGOcomparison_data %>% filter(Source == "Simulated"),
              aes(x = factor(week), y = rt),
              fill = "lightgray", color = "black", alpha = 0.6, trim = FALSE) +
  geom_jitter(data = OGOcomparison_data %>% filter(Source == "Observed"),
              aes(x = factor(week), y = rt, color = factor(subject)),
              width = 0.1, size = 1.2, alpha = 0.6) +
  facet_wrap(~ subject) +
  labs(title = "BEESTS Model Predictions vs Observed GO RTs",
       x = "Week",
       y = "Response Time (s)",
       color = "Subject") +
  coord_cartesian(ylim = c(NA, NA)) +
  theme_light()
ggsave("Images/OGOcomparison_plot1.png", OGOcomparison_plot1, width = 5, height = 5, dpi = 300)


# Repeat for STOP trials image
# Create clean subset for observed GO trials
obs_stop_subset <- OSARIobsStop %>%
  select(subject, week, rt) %>%
  mutate(Source = "Observed")

# Create clean subset for simulated GO trials
sim_stop_subset <- exgSSsimStop %>%
  select(subject, week, RT) %>%
  rename(rt = RT) %>%
  mutate(Source = "Simulated")

# Combine both into one comparison dataset
OSTOPcomparison_data <- rbind(obs_stop_subset, sim_stop_subset)

# Overlay violin plot where Simulated = violin, Observed = jitter
# Boxes split by subject, with 3 weeks inside
OSTOPcomparison_plot <- ggplot() +
  geom_violin(data = OSTOPcomparison_data %>% filter(Source == "Simulated"),
              aes(x = factor(week), y = rt),
              fill = "lightgray", color = "black", alpha = 0.6, trim = FALSE) +
  geom_jitter(data = OSTOPcomparison_data %>% filter(Source == "Observed"),
              aes(x = factor(week), y = rt, color = factor(subject)),
              width = 0.1, size = 1.2, alpha = 0.6) +
  facet_wrap(~ subject) +
  labs(title = "BEESTS Model Predictions vs Observed STOP RTs",
       x = "Week",
       y = "Response Time (s)",
       color = "Subject") +
  coord_cartesian(ylim = c(NA, NA)) +
  theme_light()
ggsave("Images/OSTOPcomparison_plot.png", OSTOPcomparison_plot, width = 5, height = 5, dpi = 300)


#Percent HDI plots
#NOTE: Only use in posters, report data instead

#Make violin plots of the HDIresults for GO Trials
OGOpctHDI_plot <- ggplot(OSARI_HDIcheck_GO, aes(x = factor(Week), y = Percent)) + #sets up the data
  geom_violin(fill = "lightgray", color = "black", alpha = 0.6)  + #sets the violin color and opacity
  geom_point(aes(color = factor(Subject))) + #adds subjects as points, with diff colors
  geom_line(aes(group = Subject, color = factor(Subject))) + #adds line across weeks for each subject
  geom_hline(yintercept = 66, linetype = "dashed", color = "red") + #adds line showing cut-off for HDI
  coord_cartesian(ylim = c(NA, 73)) + #makes it stop zooming in on the plots so much
  labs(title = "BEESTS Model GO RT Prediction Coverage",
       x = "Week",
       y = "Percent of Observed RTs Within 66% HDI",
       color = "Subject") +
  theme_light() #+  #makes it so the black background turns white

ggsave("Images/OGOpctHDI_plot.png", plot = OGOpctHDI_plot, width = 5, height = 5, dpi= 300)

#Make violin plots of the HDIresults for STOP Trials
OSTOPpctHDI_plot <- ggplot(OSARI_HDIcheck_GO, aes(x = factor(Week), y = Percent)) + #sets up the data
  geom_violin(fill = "lightgray", color = "black", alpha = 0.6)  + #sets the violin color and opacity
  geom_point(aes(color = factor(Subject))) + #adds subjects as points, with diff colors
  geom_line(aes(group = Subject, color = factor(Subject))) + #adds line across weeks for each subject
  geom_hline(yintercept = 66, linetype = "dashed", color = "red") + #adds line showing cut-off for HDI
  coord_cartesian(ylim = c(NA, 73)) + #makes it stop zooming in on the plots so much
  labs(title = "BEESTS Model GO RT Prediction Coverage",
       x = "Week",
       y = "Percent of Observed RTs Within 66% HDI",
       color = "Subject") +
  theme_light() #+  #makes it so the black background turns white

ggsave("Images/OSTOPpctHDI_plot.png", plot = OSTOPpctHDI_plot, width = 5, height = 5, dpi= 300)




### Quantile Correlations

#GO TRIALS
#Make combined data set for comparison, but need to have the source first
OGOobs_comp <- OSARIobsGo %>%
  select(subject, week, rt) %>%
  mutate(Source = "Observed")


OGOsim_comp <- exgSSsimGo %>%
  select(subject, week, RT) %>%
  rename(rt = RT) %>%
  mutate(Source = "Simulated")

OGO_comp <- bind_rows(OGOobs_comp, OGOsim_comp)


# Define quantile breaks
quantile_bins <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)

# Assign quantile bin labels 
bin_labels <- c("0-10%", "10-30%", "30-50%", "50-70%", "70-90%", "90-100%")

# Add quantile bins *within each subject-week*
OGOcomparison_binned <- OGO_comp %>%
  group_by(subject, week, Source) %>%
  mutate(quantile_bin = cut(rt, breaks = quantile(rt, probs = quantile_bins, na.rm = TRUE), 
                            include.lowest = TRUE, labels = bin_labels)) %>%
  ungroup()


# Find the match between observed and simulated
# Average RT by subject, week, source, and quantile bin
OGO_bin_summary <- OGOcomparison_binned %>%
  group_by(subject, week, quantile_bin, Source) %>%
  summarise(mean_rt = mean(rt), .groups = "drop") %>%
  pivot_wider(names_from = Source, values_from = mean_rt, names_prefix = "rt_") %>%
  filter(!is.na(rt_Observed) & !is.na(rt_Simulated))

# Run the correlations
GO_cor_results <- OGO_bin_summary %>%
  group_by(quantile_bin) %>%
  summarise(correlation = cor(rt_Observed, rt_Simulated),
            p_value = cor.test(rt_Observed, rt_Simulated)$p.value,
            n = n())

write.csv(GO_cor_results, file = "Output/OSARIGO_Hyp1_Correlations.csv", row.names = FALSE)

#STOP TRIALS
#Make combined data set for comparison, but need to have the source first
OSTOPobs_comp <- OSARIobsStop %>%
  select(subject, week, rt) %>%
  mutate(Source = "Observed")


OSTOPsim_comp <- exgSSsimStop %>%
  select(subject, week, RT) %>%
  rename(rt = RT) %>%
  mutate(Source = "Simulated")

OSTOP_comp <- bind_rows(OSTOPobs_comp, OSTOPsim_comp) %>%
  filter(!is.na(rt)) #Drop NAs, cause too small of obs to handle it


# Define quantile breaks
quantile_bins <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0) #Technically a duplicate, but keeping for consistency

# Assign quantile bin labels 
bin_labels <- c("0-10%", "10-30%", "30-50%", "50-70%", "70-90%", "90-100%")


OSTOPcomparison_binned <- OSTOP_comp %>%
  filter(!is.na(rt)) %>%
  group_by(Source) %>%
  mutate(quantile_bin = cut(rt,
                            breaks = quantile(rt, probs = quantile_bins, na.rm = TRUE),
                            include.lowest = TRUE,
                            labels = bin_labels)) %>%
  ungroup()


# Find the match between observed and simulated
# Average RT by subject, week, source, and quantile bin
OSTOP_bin_summary <- OSTOPcomparison_binned %>%
  group_by(subject, week, quantile_bin, Source) %>%
  summarise(mean_rt = mean(rt), .groups = "drop") %>%
  pivot_wider(names_from = Source, values_from = mean_rt, names_prefix = "rt_") %>%
  filter(!is.na(rt_Observed) & !is.na(rt_Simulated))

# Then do  the correlations
STOP_cor_results <- OSTOP_bin_summary %>%
  group_by(quantile_bin) %>%
  summarise(correlation = cor(rt_Observed, rt_Simulated),
            p_value = cor.test(rt_Observed, rt_Simulated)$p.value,
            n = n())

write.csv(STOP_cor_results, file = "Output/OSARISTOP_Hyp1_Correlations.csv", row.names = FALSE)


