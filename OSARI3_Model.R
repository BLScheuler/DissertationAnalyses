#OSARI Step 3
#Set up eGSS model fit

#Fitting OSARI data with exgSS model

library(tidyr)
library(dplyr)
require(msm)  # For truncated normal priors 
require(coda) # Sampling output analysis
require(loo) # For WAIC and looaic calculation
require(hypergeo) # For population plausible values
require(statmod) # Wald model
require(pracma)  # For gng and stop signal robust integration
require(numDeriv) # Prior transformations
require(vioplot) # Stop signal graphs

#set working directory

rm(list=ls())
source ("dmc/file_utils.R")
source("dmc/dmc.R")
# Load in all the DMC modules
source ("dmc/dmc_model.R")
source ("dmc/dmc_sampling.R")
source ("dmc/dmc_hierarchical.R")
source ("dmc/dmc_plotting.R")
source ("dmc/dmc_analysis.R")


load_model("EXG-SS", "exgSS.R")


all_dat <- read.csv("osari.csv")
all_dat$subject <- as.numeric(factor(all_dat$subject))



subjects <- unique(all_dat$subject)
weeks <- unique(all_dat$week)

all_dat <- subset(all_dat, rt < 2000)
all_dat$rt[all_dat$rt >= 1000] <- NA
#Changed >= to >
all_dat$rt <- all_dat$rt / 1000

model <- model.dmc(
  # SS stands for trial type (GO or Stop-signal [SS]):
  factors=list(S= "s1", SS=c("GO","SS")),
  # NR stands for "No response", i.e., go omission & successful inhibitions:
  responses=c("NR","r1"),
  # Match scores correct responses for each GO stimulus as usual, and scores 
  # the "correct" stimulus corresponding to an NR response, but the latter has
  # no effect (it just avoids a standard check making sure that each response is scored): 
  match.map=list(M=list(s1="r1",s1="NR")), #same s1 twice
  p.map=list(mu="M",sigma="M",tau="M",tf="1",muS="1",sigmaS="1",tauS="1",gf="1"),
  # No errors:
  constants=c(mu.false=1e6,sigma.false=.001,tau.false=.001, tf=0, gf=0),
  type="exgss")

# This gives mean GoRT of .5 + .08 = .58 and SD Go RT of sqrt(.05^2+.08^2) = 0.09 
# and mean SSRT of .2+.05 = 0.25 and SD SSRT of sqrt(.03^2+.05^2) = 0.06:
p.vector  <- c(mu.true=.5,muS=.2,sigma.true=.05,sigmaS=.03,tau.true=.08,tauS=.05) 
# Uniform (scaled beta) priors:
p1 <- p.vector; p1[1:length(p1)] <- 1; p1
p.prior <- prior.p.dmc(
  dists = rep("beta",length(p1)),p1=p1,p2=rep(1,length(p1)), # Uniform(0,1)
  lower=rep(0,length(p1)),upper=c(1,1,rep(.5,4)) # Scale to Uniform(lower,upper)
)


samples_week <- vector("list", length(weeks))
samples_subject_week <- vector("list", length(subjects))
for (sj in subjects) { 
  samples_subject_week[[sj]] <- samples_week
}


for (sj in subjects) {
  for (wk in weeks) {
    #dat <- all_dat[all_dat$subject == sj & all_dat$week == wk,]
    dat <- subset(all_dat, subject == sj & week == wk)
    
    is_outlier <- with(dat, (rt > mean(rt, na.rm=TRUE) + 3 * sd(rt, na.rm=TRUE)) |  (rt < mean(rt, na.rm=TRUE) - 3 * sd(rt, na.rm=TRUE)))
    dat <- subset(dat, !is_outlier)
    
    dat_dmc <- data.frame(S  = "s1",
                          SS = ifelse(dat$ss_presented, "SS", "GO"), 
                          R  = ifelse(dat$inhibited == 1, "NR", "r1"),
                          RT = dat$rt,
                          SSD= ifelse(dat$ssd == .500, .5, Inf))

    
    data_dmc <- data.model.dmc(dat_dmc,model)
    
    # Start sampling
    samples <- samples.dmc(nmc=100,p.prior,data_dmc)
    samples <- run.unstuck.dmc(samples,report = 10,cores=4,p.migrate=.05,verbose=TRUE)
    
    samples2 <- run.converge.dmc(samples.dmc(samples=samples,nmc=50,thin=5),
                                 report=10, cores=4,cut=1.1,verbose=TRUE,nmc=50,minN=500,max.try=20)
    
    samples_subject_week[[sj]][[wk]] <- samples2
    save(samples_subject_week, file="osari_posterior.RData")
  }
}



#Create SSRTs (originally from plotOSARIparms script)
load("osari_posterior.RData")

group_by_week <- function(dat, parameters) { 
  weeks <- 1:length(dat[[1]])
  subjects <- 1:length(dat)
  samp_dim <- dim(dat[[1]][[1]]$theta[,1,])
  out_df <- matrix(NA, 0, 3 + length(parameters))
  colnames(out_df) <- c("subject", "week", "sample", parameters)
  for (wk in weeks) { 
    for (sj in subjects) { 
      out_mat <- matrix(NA, prod(samp_dim), 0)
      for (parameter in parameters) { 
        out_mat <- cbind(out_mat, c(dat[[sj]][[wk]]$theta[,parameter,1:samp_dim[2]]))
      }
      out_df <- rbind(out_df, cbind(sj, wk, 1:prod(samp_dim), out_mat))
    }
  }
  out_df <- data.frame(out_df)
  out_df$subject <- factor(out_df$subject)
  out_df$week <- factor(out_df$week, ordered=TRUE)
  return(out_df)
}

parameters <- c("muS", "tauS","sigmaS", "mu.true", "tau.true", "sigma.true")
beests_df <- group_by_week(samples_subject_week, parameters)

beests_subj_means <- beests_df %>% 
  group_by(subject, week) %>% 
  summarise(muS = mean(muS), tauS=mean(tauS), sigmaS = mean(sigmaS), mu.true=mean(mu.true), tau.true=mean(tau.true), sigma.true=mean(sigma.true))

beests_week_means <- beests_df %>% 
  group_by(week, sample) %>% 
  summarise(muS = mean(muS), tauS=mean(tauS), sigmaS = mean(sigmaS), mu.true=mean(mu.true), tau.true=mean(tau.true), sigma.true=mean(sigma.true))

#Combine to make SSRT for week group means
beests_week_means$SSRT <- with(beests_week_means, muS + tauS)

#Combine mu and tau to make SSRT for individuals
beests_subj_means$SSRT <- with(beests_subj_means, muS + tauS)

save(beests_subj_means, file = "BeestsSubjMeans.RData")


