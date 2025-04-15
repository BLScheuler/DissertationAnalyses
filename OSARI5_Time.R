#OSARI Step 5
#Change over time


### SET UP THE MODEL

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)

#set working directory

#load exgSS modeling data
load("BeestsSubjMeans.RData")


#Set up data for ulam model
ECmodData <- list(
  part = as.integer(as.factor(beests_subj_means$subject)),
  Week2 = 1*(beests_subj_means$week == 2),
  Week3 = 1*(beests_subj_means$week == 3),
  SSRT = beests_subj_means$SSRT
)


#Set seed for reproducibility
set.seed(123)

#Model including only time
m_TimeSSRT <- ulam(
  alist(
    ## Time -> Exec. Con. 
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
    Rho_SSRT ~ lkj_corr(2)
    
  ), data = ECmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_TimeSSRT, depth=3)


### EVALUATE OUTPUT FOR CHANGE OVER TIME

#Evaluate if the model posteriors are in the predicted directions
#Extract samples from the model
post <- extract.samples(m_TimeSSRT)


#Change from wk 1 to wk2 (b1_SSRTwk is the slope between wk1 and wk2)
wk1towk2 <- mean(post$b1_SSRTwk)
#-0.01712206

#Proportion of posterior samples where SSRT increased from wk1 to wk2
mean(post$b1_SSRTwk > 0)
#0.0985

#View where the middle 66% of the posterior distribution lies
HPDI(post$b1_SSRTwk, prob = 0.66)
#-0.02823850 -0.00609111 

#Change from wk 2 to wk3 (b2_SSRTwk is the slope between wk2 and wk3)
wk2towk3 <- mean(post$b2_SSRTwk)
#-0.008333134

#Proportion of posterior samples where SSRT decreased from wk2 to wk3
mean(post$b2_SSRTwk < 0)
#0.6814

#View where the middle 66% of the posterior distribution lies
HPDI(post$b2_SSRTwk, prob = 0.66)
#-0.0274470  0.0090224 

#Save that all in a database for later
SSRT_summary <- data.frame(
  Weeks = c("Wk1toWk2", "Wk2toWk3"),
  Slope = c(mean(post$b1_SSRTwk), mean(post$b2_SSRTwk)),
  Proportion = c(mean(post$b1_SSRTwk > 0), mean(post$b2_SSRTwk < 0)),  
  LowHDI = c(HPDI(post$b1_SSRTwk, prob = 0.66)[1], HPDI(post$b2_SSRTwk, prob = 0.66)[1]),
  HighHDI = c(HPDI(post$b1_SSRTwk, prob = 0.66)[2], HPDI(post$b2_SSRTwk, prob = 0.66)[2])
)

# Save to CSV
write.csv(SSRT_summary, "Output/Hyp2_SSRT_Summary.csv", row.names = FALSE)




### PLOT FOR VISUAL

# --- HDIs ---
hdi1 <- HPDI(post$b1_SSRTwk, prob = 0.66)
hdi2 <- HPDI(post$b2_SSRTwk, prob = 0.66)

# --- Format Posterior Samples for Plotting ---
post_df <- data.frame(
  Slope = c(post$b1_SSRTwk, post$b2_SSRTwk),
  Comparison = factor(
    rep(c("Week 1 to Week 2", "Week 2 to Week 3"), each = length(post$b1_SSRTwk)),
    levels = c("Week 1 to Week 2", "Week 2 to Week 3")
  )
)

# --- HDI Segment Data ---
dens_wk1wk2 <- density(post$b1_SSRTwk)
dens_wk2wk3 <- density(post$b2_SSRTwk)


hdi_df <- data.frame(
  Comparison = c("Week 1 to Week 2", "Week 2 to Week 3"),
  x = c(hdi1[1], hdi2[1]),
  xend = c(hdi1[2], hdi2[2]),
  y = c(max(dens_wk1wk2$y), max(dens_wk2wk3$y)) * 0.3  # places HDI line about a third of the way up
)

# --- Custom Legend Guide Data ---
legend_df <- data.frame(
  x = c(NA, NA),
  xend = c(NA, NA),
  y = c(NA, NA),
  yend = c(NA, NA),
  type = c("Null Effect", "66% HDI")
)

# --- Plot ---
SSRT_timePlot <- ggplot(post_df, aes(x = Slope)) +
  geom_density(aes(fill = Comparison), alpha = 0.6, show.legend = FALSE) +
  
  # Dashed line for null effect (manually add to legend)
  geom_vline(aes(xintercept = 0, linetype = "Null Effect"), color = "black", linewidth = 1) +
  
  # HDI bar (manually add to legend)
  geom_segment(
    data = hdi_df,
    aes(x = x, xend = xend, y = y, yend = y, linetype = "66% HDI"),
    color = "black",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(~Comparison, scales = "free") +
  
  scale_linetype_manual(
    name = "Legend", 
    values = c("Null Effect" = "dashed", "66% HDI" = "solid")
  ) +
  
  labs(
    title = "Posterior Distributions of SSRT Slope",
    x = "Posterior Slope (SSRT Change)",
    y = "Density"
  ) +
  
  theme_light() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.box.background = element_rect(color = "black", size = 0.8),
    legend.box.margin = margin(5, 5, 5, 5),
    legend.title = element_blank(), #Hides "Legend"
    legend.key.height = unit(1.8, "lines") #makes it so we can see the dashed line in saved image
  )

#Save the image
ggsave("Images/Hyp2_SSRTtimePlot.png", SSRT_timePlot, width = 8, height = 5, dpi = 300)





