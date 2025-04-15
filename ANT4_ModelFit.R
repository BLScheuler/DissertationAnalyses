#ANT Step 4
#Test model fit

#Set working directory


set.seed(123)


library(coda)     # For computing HPD intervals (HDI)
library(dplyr)    # For data manipulation
library(ggplot2)  # For the violin plots
library(flankr)   # Load packages specific to SSP
library(tidyr)    # For the correlations
library(purrr)    # For the correlations


# Load ANT data and best-fitting model parameters
ANTdata <- read.csv("CleanANT.csv")
ANTpar_bf <- read.csv("ANTparameters.csv")


# Function to simulate RTs using Shrinking Spotlight Model
simulate_rt_SSP <- function(A, Ter, p, rd, sda, n = 10000) {
  simulated <- simulateSSP(parms = c(A, Ter, p, rd, sda), nTrials = n)
  return(simulated$rt)
}

# Data frame to store HDI results
ANT_HDIcheck <- data.frame(subject = numeric(),
                           week = numeric(),
                           pct_in_HDI = numeric())


## Check overall coverage of RTs
# Loop through each subject-week combo
for (i in 1:nrow(ANTpar_bf)) {
  
  subj_id <- ANTpar_bf$subject[i]
  curr_week <- ANTpar_bf$week[i]
  
  # Actual observed RTs from the ANT data
  ANTobserved_rts <- subset(ANTdata, subject == subj_id & week == curr_week)$rt
  
  # Simulated RTs from the best-fitting model
  sim_rts <- simulate_rt_SSP(A   = ANTpar_bf$A[i],
                             Ter = ANTpar_bf$Ter[i],
                             p   = ANTpar_bf$p[i],
                             rd  = ANTpar_bf$rd[i],
                             sda = ANTpar_bf$sda[i])
  
  # Get 66% HDI of simulated RTs
  mcmc_sim <- as.mcmc(sim_rts)  #coda package needs RTs as MCMC object before it can do HPD interval
  ANT_HDIbounds <- HPDinterval(mcmc_sim, prob = 0.66)
  lower <- ANT_HDIbounds[1]
  upper <- ANT_HDIbounds[2]
  
  # What percent of observed RTs fall within the HDI
  in_hdi <- ANTobserved_rts >= lower & ANTobserved_rts <= upper
  pct <- mean(in_hdi) * 100
  
  # Save results
  ANT_HDIcheck[nrow(ANT_HDIcheck) + 1, ] <- c(subj_id, curr_week, pct)
}

# Save results to CSV
write.csv(ANT_HDIcheck, file = "Output/ANT_Hyp1_Overallresults.csv", row.names = FALSE)


####

#Now do quantile specific HDI analyses

#Set up quantile materials
quantile_levels <- c(0.1, 0.3, 0.5, 0.7, 0.9)
n_simulations <- 2000
HDI_QuantileResults <- data.frame()

# Loop through each subject-week combo
for (row in 1:nrow(ANTpar_bf)) {
  subj_id <- ANTpar_bf$subject[row]
  week_id <- ANTpar_bf$week[row]
  
  #Pull in observed RTs
  obs_rts <- subset(ANTdata, subject == subj_id & week == week_id)$rt
  n_obs <- length(obs_rts)
  if (n_obs < 10) next
  
  #Create play to store quantiles
  sim_quantiles <- matrix(NA, nrow = n_simulations, ncol = length(quantile_levels))
  
  #Simulations loop for parameters
  for (i in 1:n_simulations) {
    sim_rts <- simulate_rt_SSP(
      A = ANTpar_bf$A[row],
      Ter = ANTpar_bf$Ter[row],
      p = ANTpar_bf$p[row],
      rd = ANTpar_bf$rd[row],
      sda = ANTpar_bf$sda[row],
      n = n_obs
    )
    if (all(is.finite(sim_rts))) {
      sim_quantiles[i, ] <- quantile(sim_rts, probs = quantile_levels)
    }
  }
  
  #Simulation quantiles
  sim_quantiles <- sim_quantiles[complete.cases(sim_quantiles), ]
  if (nrow(sim_quantiles) < 50) next
  
  #Observed quantiles
  obs_quantiles <- quantile(obs_rts, probs = quantile_levels)
  HDIs <- apply(sim_quantiles, 2, function(qvals) HPDinterval(as.mcmc(qvals), prob = 0.66))
  
  #check if they are within the HDI
  within_HDI <- mapply(function(obs, bounds) obs >= bounds[1] & obs <= bounds[2],
                       obs = obs_quantiles,
                       bounds = as.data.frame(HDIs))
  
  # Store full results per quantile
  result_row <- data.frame(
    subject = subj_id,
    week = week_id,
    q10_obs = obs_quantiles[1], # obs rt at q10
    q10_HDI_lower = HDIs[1,1],  # lower HDI from sim for q10
    q10_HDI_upper = HDIs[2,1],  # upper HDI from sim for q10
    q10_in_HDI = within_HDI[1], # whether or not obs was in HDI
    q30_obs = obs_quantiles[2],
    q30_HDI_lower = HDIs[1,2],
    q30_HDI_upper = HDIs[2,2],
    q30_in_HDI = within_HDI[2],
    q50_obs = obs_quantiles[3],
    q50_HDI_lower = HDIs[1,3],
    q50_HDI_upper = HDIs[2,3],
    q50_in_HDI = within_HDI[3],
    q70_obs = obs_quantiles[4],
    q70_HDI_lower = HDIs[1,4],
    q70_HDI_upper = HDIs[2,4],
    q70_in_HDI = within_HDI[4],
    q90_obs = obs_quantiles[5],
    q90_HDI_lower = HDIs[1,5],
    q90_HDI_upper = HDIs[2,5],
    q90_in_HDI = within_HDI[5],
    num_in_HDI = sum(within_HDI)
  )
  
  #Stores quantile results by adding each row (1 subject + 1 week) below each other
  HDI_QuantileResults <- rbind(HDI_QuantileResults, result_row)
}
write.csv(HDI_QuantileResults, file = "Output/ANT_Hyp1_Quantileresults.csv", row.names = FALSE)



### PLOTS

# See what the actual RTs looked like
# First, calculate subject-level mean RTs per week (for lines)
ANTdata_summary <- ANTdata %>%
  group_by(subject, week) %>%
  summarise(mean_rt = mean(rt), .groups = "drop")

# Create the plot
ANTrt_plot <- ggplot(ANTdata, aes(x = factor(week), y = rt)) +
  geom_violin(fill = "lightgray", color = "black", alpha = 0.6, trim = FALSE) +
  geom_jitter(aes(color = factor(subject)), width = 0.1, size = 1.2, alpha = 0.6) +
  labs(title = "ANT Observed RT Distributions",
       x = "Week",
       y = "Observed Response Time (s)",
       color = "Subject") +
  coord_cartesian(ylim = c(0, NA)) +  # allows auto-scaling for upper bound
  theme_light()


# Save the file 
ggsave("Images/ANTrt_plot.png", ANTrt_plot, width = 5, height = 5, dpi = 300)


#Create a plot overlaying the model on the observed RTs

# Create an empty data frame to store all simulated RTs
ANTsimulated_data <- data.frame()

# Loop through each subject-week to simulate RTs using best-fitting parameters
for (i in 1:nrow(ANTpar_bf)) {
  subj_id <- ANTpar_bf$subject[i]
  curr_week <- ANTpar_bf$week[i]
  
  # Simulate 1000 RTs for that subject-week
  sim_rts <- simulate_rt_SSP(A = ANTpar_bf$A[i],
                             Ter = ANTpar_bf$Ter[i],
                             p = ANTpar_bf$p[i],
                             rd = ANTpar_bf$rd[i],
                             sda = ANTpar_bf$sda[i],
                             n = 1000)
  
  # Store the results with labels
  temp_df <- data.frame(
    subject = subj_id,
    week = curr_week,
    rt = sim_rts,
    source = "Simulated"
  )
  
  ANTsimulated_data <- rbind(ANTsimulated_data, temp_df)
}

# Prepare observed data
ANTobserved_data <- ANTdata %>%
  dplyr::select(subject, week, rt) %>%
  mutate(source = "Observed")

# Combine both into one dataset
ANTcomparison_data <- rbind(ANTobserved_data, ANTsimulated_data)
ANTcomparison_data$source <- factor(ANTcomparison_data$source, levels = c("Observed", "Simulated"))

# Overlay violin plot
ANTcomparison_plot <- ggplot() +
  geom_violin(data = ANTcomparison_data %>% filter(source == "Simulated"),
              aes(x = factor(week), y = rt),
              fill = "lightgray", color = "black", alpha = 0.6, trim = FALSE) +
  geom_jitter(data = ANTcomparison_data %>% filter(source == "Observed"),
              aes(x = factor(week), y = rt, color = factor(subject)),
              width = 0.1, size = 1.2, alpha = 0.6) +
  facet_wrap(~ subject) +
  labs(title = "Shrinking Spotlight Model Predictions vs Observed RTs",
       x = "Week",
       y = "Response Time (s)",
       color = "Subject") +
  coord_cartesian(ylim = c(0, NA)) +
  theme_light()

ggsave("Images/ANTcomparison_plot.png", ANTcomparison_plot, width = 5, height = 5, dpi = 300)



#Percent HDI plots
#NOTE: Only use in posters, report data instead

#Make violin plots of the HDIresults
ANTpctHDI_plot <- ggplot(ANT_HDIcheck, aes(x = factor(week), y = pct_in_HDI)) + #sets up the data
  geom_violin(fill = "lightgray", color = "black", alpha = 0.6)  + #sets the violin color and opacity
  geom_point(aes(color = factor(subject))) + #adds subjects as points, with diff colors
  geom_line(aes(group = subject, color = factor(subject))) + #adds line across weeks for each subject
  geom_hline(yintercept = 66, linetype = "dashed", color = "red") + #adds line showing cut-off for HDI
  coord_cartesian(ylim = c(60, 80)) + #makes it stop zooming in on the plots so much
  labs(title = "Shrinking Spotlight Model RT Prediction Coverage",
       x = "Week",
       y = "Percent of Observed RTs Within 66% HDI",
       color = "Subject") +
  coord_cartesian(ylim = c(60, NA)) +
  theme_light() #+  #makes it so the black background turns white

ggsave("Images/ANTpctHDI_plot.png", plot = ANTpctHDI_plot, width = 5, height = 5, dpi= 300)



### Quantile Correlations

# Define quantile breaks
quantile_bins <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0)

# Assign quantile bin labels 
bin_labels <- c("0-10%", "10-30%", "30-50%", "50-70%", "70-90%", "90-100%")

# Add quantile bins *within each subject-week*
ANTcomparison_binned <- ANTcomparison_data %>%
  group_by(subject, week, source) %>%
  mutate(quantile_bin = cut(rt, breaks = quantile(rt, probs = quantile_bins), include.lowest = TRUE, labels = bin_labels)) %>%
  ungroup()

# Find the match between observed and simulated
# Average RT by subject, week, source, and quantile bin
bin_summary <- ANTcomparison_binned %>%
  group_by(subject, week, quantile_bin, source) %>%
  summarise(mean_rt = mean(rt), .groups = "drop") %>%
  pivot_wider(names_from = source, values_from = mean_rt, names_prefix = "rt_") %>%
  filter(!is.na(rt_Observed) & !is.na(rt_Simulated))

# Run the correlations
cor_results <- bin_summary %>%
  group_by(quantile_bin) %>%
  summarise(correlation = cor(rt_Observed, rt_Simulated),
            p_value = cor.test(rt_Observed, rt_Simulated)$p.value,
            n = n())

write.csv(cor_results, file = "Output/ANT_Hyp1_Correlations.csv", row.names = FALSE)


