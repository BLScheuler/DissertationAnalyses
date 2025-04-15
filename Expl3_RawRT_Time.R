# Exploratory Step 3
#Raw RT over time


### SET UP THE MODEL

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)

#set working directory

#Pull in OSARI raw data
OSARIdata <- read.csv("osari.csv")

#Keep only RTs that are from unsuccessful stop trials
OSARIdata <- OSARIdata[OSARIdata$ss_presented == 1 & OSARIdata$inhibited == 0, ]


#Get the average rawS for each participant and week
MeanRTdata <- OSARIdata %>%
  group_by(subject, week) %>%
  summarize(meanRT = mean(rt, na.rm = TRUE), .groups = "drop")

# Convert from ms to s
MeanRTdata$meanRT <- MeanRTdata$meanRT / 1000


#Set up data for ulam model
ECmodData <- list(
  part = as.integer(as.factor(MeanRTdata$subject)),
  Week2 = 1*(MeanRTdata$week == 2),
  Week3 = 1*(MeanRTdata$week == 3),
  RT = MeanRTdata$meanRT
)


#Set seed for reproducibility
set.seed(123)

#Model including only time
m_TimeRT <- ulam(
  alist(
    ## Time -> Exec. Con. 
    #distribution for EC parameter 
    RT ~ dnorm(mu_RT, sigma_RT),
    #Set up participant change over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3
    mu_RT <- p_RT[part] + p_RTwk2[part]*Week2 + p_RTwk3[part]*Week3,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_RT,p_RTwk2, p_RTwk3)[part] ~ multi_normal( c(a_RT,b1_RTwk, b2_RTwk) , Rho_RT , sigma_indRT),
    
    #RT priors
    a_RT ~ normal(.8,.2),
    b1_RTwk ~ normal(0,.1),
    b2_RTwk ~ normal(0,.1),
    
    
    sigma_RT ~ exponential(1),
    sigma_indRT ~ exponential(1),
    Rho_RT ~ lkj_corr(2)
    
  ), data = ECmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_TimeRT, depth=3)


### EVALUATE OUTPUT FOR CHANGE OVER TIME

#Evaluate if the model posteriors are in the predicted directions
#Extract samples from the model
post <- extract.samples(m_TimeRT)


#Change from wk 1 to wk2 (b1_RTwk is the slope between wk1 and wk2)
wk1towk2 <- mean(post$b1_RTwk)
#-0.03256809

#Proportion of posterior samples where RT increased from wk1 to wk2
mean(post$b1_RTwk > 0)
#0.0974

#View where the middle 66% of the posterior distribution lies
HPDI(post$b1_RTwk, prob = 0.66)
#-0.02823850 -0.00609111 

#Change from wk 2 to wk3 (b2_RTwk is the slope between wk2 and wk3)
wk2towk3 <- mean(post$b2_RTwk)
#-0.0593811 -0.0134487

#Proportion of posterior samples where RT decreased from wk2 to wk3
mean(post$b2_RTwk < 0)
# 0.8034

#View where the middle 66% of the posterior distribution lies
HPDI(post$b2_RTwk, prob = 0.66)
#-0.04090620  0.00249431 

#Save that all in a database for later
RT_summary <- data.frame(
  Weeks = c("Wk1toWk2", "Wk2toWk3"),
  Slope = c(mean(post$b1_RTwk), mean(post$b2_RTwk)),
  Proportion = c(mean(post$b1_RTwk > 0), mean(post$b2_RTwk < 0)),  
  LowHDI = c(HPDI(post$b1_RTwk, prob = 0.66)[1], HPDI(post$b2_RTwk, prob = 0.66)[1]),
  HighHDI = c(HPDI(post$b1_RTwk, prob = 0.66)[2], HPDI(post$b2_RTwk, prob = 0.66)[2])
)

# Save to CSV
write.csv(RT_summary, "Output/Expl_RT_TimeSumm.csv", row.names = FALSE)




### PLOT FOR VISUAL

# --- HDIs ---
hdi1 <- HPDI(post$b1_RTwk, prob = 0.66)
hdi2 <- HPDI(post$b2_RTwk, prob = 0.66)

# --- Format Posterior Samples for Plotting ---
post_df <- data.frame(
  Slope = c(post$b1_RTwk, post$b2_RTwk),
  Comparison = factor(
    rep(c("Week 1 to Week 2", "Week 2 to Week 3"), each = length(post$b1_RTwk)),
    levels = c("Week 1 to Week 2", "Week 2 to Week 3")
  )
)

# --- HDI Segment Data ---
dens_wk1wk2 <- density(post$b1_RTwk)
dens_wk2wk3 <- density(post$b2_RTwk)


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
RT_timePlot <- ggplot(post_df, aes(x = Slope)) +
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
    title = "Posterior Distributions of RT Slope",
    x = "Posterior Slope (RT Change)",
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
ggsave("Images/Expl_RTtimePlot.png", RT_timePlot, width = 8, height = 5, dpi = 300)
