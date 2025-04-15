#DNB Step 5
#Change over time assessment

### SET UP THE MODEL

#load necessary packages
library(rethinking)
library(ggplot2)
library(coda)

#set working directory

#load DNB cap coef data data
DNBdata <- read.csv("DNB_CapCoef.csv")

#Set up data for ulam model
WCmodData <- list(
  part = as.integer(as.factor(DNBdata$Subject)),
  Week2 = 1*(DNBdata$Week == 2),
  Week3 = 1*(DNBdata$Week == 3),
  WC = DNBdata$Cz
)

#Set seed for reproducibility
set.seed(123)

#Model including only time
m_TimeWC <- ulam(
  alist(
    ## Time -> Workload Capacity 
    #distribution for WC parameter 
    WC ~ dnorm(mu_WC, sigma_WC),
    #Set up participant change over time
    #Week 1 (where wk2 and wk3 ar zero) plus changes for wk 2 and wk 3
    mu_WC <- p_WC[part] + p_WCwk2[part]*Week2 + p_WCwk3[part]*Week3,
    
    #multivariate normal priors where individuals can experience different changes per week
    c(p_WC,p_WCwk2, p_WCwk3)[part] ~ multi_normal( c(a_WC,b1_WCwk, b2_WCwk) , Rho_WC , sigma_indWC),
    
    #WC Priors
    a_WC ~ normal(3,.75),
    b1_WCwk ~ normal(0,.5),
    b2_WCwk ~ normal(0,.5),
    
    
    sigma_WC ~ exponential(1),
    sigma_indWC ~ exponential(1),
    Rho_WC ~ lkj_corr(2)
    
  ), data = WCmodData, chains = 4, iter = 5000, log_lik = TRUE)

precis(m_TimeWC, depth=3)

### EVALUATE OUTPUT FOR CHANGE OVER TIME

#Evaluate if the model posteriors are in the predicted directions
#Extract samples from the model
post <- extract.samples(m_TimeWC)


#Change from wk 1 to wk2 (b1_WCwk is the slope between wk1 and wk2)
wk1towk2 <- mean(post$b1_WCwk)
#0.4752089

#Proportion of posterior samples where WC decreased from wk1 to wk2
mean(post$b1_WCwk < 0)
#0.1072

#View where the middle 66% of the posterior distribution lies
HPDI(post$b1_WCwk, prob = 0.66)
#0.185713 0.904581 

#Change from wk 2 to wk3 (b2_WCwk is the slope between wk2 and wk3)
wk2towk3 <- mean(post$b2_WCwk)
#0.1683107

#Proportion of posterior samples where WC increased from wk2 to wk3
mean(post$b2_WCwk > 0)
#0.6363

#View where the middle 66% of the posterior distribution lies
HPDI(post$b2_WCwk, prob = 0.66)
#-0.133471  0.541264 

#Save that all in a database for later
WC_summary <- data.frame(
  Weeks = c("Wk1toWk2", "Wk2toWk3"),
  Slope = c(mean(post$b1_WCwk), mean(post$b2_WCwk)),
  Proportion = c(mean(post$b1_WCwk < 0), mean(post$b2_WCwk > 0)),  
  LowHDI = c(HPDI(post$b1_WCwk, prob = 0.66)[1], HPDI(post$b2_WCwk, prob = 0.66)[1]),
  HighHDI = c(HPDI(post$b1_WCwk, prob = 0.66)[2], HPDI(post$b2_WCwk, prob = 0.66)[2])
)

# Save to CSV
write.csv(WC_summary, "Output/Hyp2_WC_Summary.csv", row.names = FALSE)




### PLOT FOR VISUAL

# --- HDIs ---
hdi1 <- HPDI(post$b1_WCwk, prob = 0.66)
hdi2 <- HPDI(post$b2_WCwk, prob = 0.66)

# --- Format Posterior Samples for Plotting ---
post_df <- data.frame(
  Slope = c(post$b1_WCwk, post$b2_WCwk),
  Comparison = factor(
    rep(c("Week 1 to Week 2", "Week 2 to Week 3"), each = length(post$b1_WCwk)),
    levels = c("Week 1 to Week 2", "Week 2 to Week 3")
  )
)

# --- HDI Segment Data ---
# Get density peak values per comparison
dens_wk1wk2 <- density(post$b1_WCwk)
dens_wk2wk3 <- density(post$b2_WCwk)

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
WC_timePlot <- ggplot(post_df, aes(x = Slope)) +
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
    title = "Posterior Distributions of Workload Capacity Slope",
    x = "Posterior Slope (WC Change)",
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
ggsave("Images/Hyp2_WCtimePlot.png", WC_timePlot, width = 8, height = 5, dpi = 300)
