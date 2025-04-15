#ANT Step 7
#Chech time analysis without poorly fit data point


### SET UP THE MODEL

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)

#set working directory

#load DNB cap coef data data
ANTdata <- read.csv("ANTparameters.csv")

#Filter out subject 5, week 1
# Convert subject to factor first
ANTdata$subject <- as.factor(ANTdata$subject)

# Overwrite ANTdata by filtering out week 1 for subject 5
ANTdata <- ANTdata[!(ANTdata$subject == "5" & ANTdata$week == 1), ]


#Set up data for ulam model
ITmodData <- list(
  part = as.integer(as.factor(ANTdata$subject)),
  Week2 = 1*(ANTdata$week == 2),
  Week3 = 1*(ANTdata$week == 3),
  IT = standardize((ANTdata$sda)/(ANTdata$rd))
)

#Set seed for reproducibility
set.seed(123)

#Model including only time
m_TimeIT <- ulam(
  alist(
    ## TIME -> ATTENTION 
    #distribution for IT parameter 
    IT ~ dnorm(mu_IT, sigma_IT),
    #Set up participant change over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3
    mu_IT <- p_IT[part] + p_ITwk2[part]*Week2 + p_ITwk3[part]*Week3,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_IT,p_ITwk2, p_ITwk3)[part] ~ multi_normal( c(a_IT,b1_ITwk, b2_ITwk) , Rho_IT, sigma_indIT),
    
    #IT Priors
    a_IT ~ normal(0, 1),
    b1_ITwk ~ normal(0, 1),
    b2_ITwk ~ normal(0, 1),
    
    
    sigma_IT ~ exponential(1),
    sigma_indIT ~ exponential(1),
    Rho_IT ~ lkj_corr(2)
    
  ), data = ITmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_TimeIT, depth=3)



### EVALUATE OUTPUT FOR CHANGE OVER TIME

#Evaluate if the model posteriors are in the predicted directions
#Extract samples from the model
post <- extract.samples(m_TimeIT)


#Change from wk 1 to wk2 (b1_ITwk is the slope between wk1 and wk2)
wk1towk2 <- mean(post$b1_ITwk)
#-0.6088299 to -0.7867593

#Proportion of posterior samples where IT increased from wk1 to wk2
mean(post$b1_ITwk > 0)
#from 0.1178 to 0.0699


#View where the middle 66% of the posterior distribution lies
HPDI(post$b1_ITwk, prob = 0.66)
#-1.086130 -0.135241  to  -1.30915 -0.37358 

#Change from wk 2 to wk3 (b2_ITwk is the slope between wk2 and wk3)
wk2towk3 <- mean(post$b2_ITwk)
#-0.6031057 to -0.7862796

#Proportion of posterior samples where IT decreased from wk2 to wk3
mean(post$b2_ITwk < 0)
#0.8735 to 0.9322

#View where the middle 66% of the posterior distribution lies
HPDI(post$b2_ITwk, prob = 0.66)
#-1.284160 -0.298746 to -1.309360 -0.385675 

#Save that all in a database for later
IT_summary <- data.frame(
  Weeks = c("Wk1toWk2", "Wk2toWk3"),
  Slope = c(mean(post$b1_ITwk), mean(post$b2_ITwk)),
  Proportion = c(mean(post$b1_ITwk > 0), mean(post$b2_ITwk < 0)),  
  LowHDI = c(HPDI(post$b1_ITwk, prob = 0.66)[1], HPDI(post$b2_ITwk, prob = 0.66)[1]),
  HighHDI = c(HPDI(post$b1_ITwk, prob = 0.66)[2], HPDI(post$b2_ITwk, prob = 0.66)[2])
)

# Save to CSV
write.csv(IT_summary, "Output/Hyp2Check_Summary.csv", row.names = FALSE)




### PLOT FOR VISUAL

# --- HDIs ---
hdi1 <- HPDI(post$b1_ITwk, prob = 0.66)
hdi2 <- HPDI(post$b2_ITwk, prob = 0.66)

# --- Format Posterior Samples for Plotting ---
post_df <- data.frame(
  Slope = c(post$b1_ITwk, post$b2_ITwk),
  Comparison = factor(
    rep(c("Week 1 to Week 2", "Week 2 to Week 3"), each = length(post$b1_ITwk)),
    levels = c("Week 1 to Week 2", "Week 2 to Week 3")
  )
)

# --- HDI Segment Data ---
dens_wk1wk2 <- density(post$b1_ITwk)
dens_wk2wk3 <- density(post$b2_ITwk)


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
IT_timePlot <- ggplot(post_df, aes(x = Slope)) +
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
    title = "Posterior Distributions of Interference Time Slope",
    x = "Posterior Slope (IT Change)",
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
ggsave("Images/Hyp2Check_ITtimePlot.png", IT_timePlot, width = 8, height = 5, dpi = 300)
