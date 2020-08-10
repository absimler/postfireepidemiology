## Code for "Wildfire alters the disturbance impacts of an emerging forest disease via changes to host occurrence and demographic structure" ###

## Part II. Does long-term fire history influence pathogen invasion, disease intensification within stands, and stand-level mortality? #####

library(rstan)
library(rstanarm)
library(bayesplot)
library(loo)
library(tidyverse)

setwd("~/Dropbox/Data Analysis/firehistorydisease")

## load all plots (n=280) for analyses based on 2006-2007 data:
Plot = read.csv("datasets/all280Plots_20062007.csv", header=T, na.strings=c("", " "))

## Load more detailed long term plot data (90 plots) for analyses based on 2013-2014 data
LTPlot <- read.csv("datasets/longtermPlots_201014.csv", header=T, na.strings=c("", " "))



# Format/scale variables:
Plot$Pr2006 <- ifelse(Plot$totinf06>0, 1, 0) ## Identify if Pr is present at plot level
Plot$ForestAllianceType <- group_indices(Plot, ForestAllianceType)
Plot$plotindex <- group_indices(Plot, BigSurPlot)
Plot$tmax.c <- scale(Plot$tmax)
Plot$totalfire.c <- scale(Plot$totalfire)
Plot$burned1950.c <- scale(Plot$burned1950)
Plot$forest.c <- scale(Plot$ForestAllianceType)

LTPlot$Pr2010 <- ifelse(LTPlot$totinf10>0, 1, 0) # Identify if Pr is present at plot level in 2010-11
LTPlot$Pr2013 <- ifelse(LTPlot$totinf13>0, 1, 0) ## Identify if Pr is present at plot level in 2013-14
LTPlot <- subset(LTPlot, tothost10>0 & Pr2010>0) ## Subset longterm plots to only infested plots
LTPlot$tmax.c <- scale(LTPlot$tmax)
LTPlot$plotindex <- group_indices(LTPlot, BigSurPlot)
LTPlot$ForestAllianceType <- group_indices(LTPlot, ForestAllianceType)
LTPlot$totalfire.c <- scale(LTPlot$totalfire)
LTPlot$burned1950.c <- scale(LTPlot$burned1950.08)
LTPlot$forest.c <- scale(LTPlot$ForestAllianceType)


## Create simulation data for plotting marginal effects:
fire.seq <- rep(c(rep(0, 20), rep(1, 20), rep(2, 20), rep(3, 20), rep(4, 20)), 2)
b50.seq <- rep(c(-1.1035656, 0.9029173), 100)
for.seq <- c(rep(-0.8272208, 100), rep(1.2045496, 100))
nsims <- length(fire.seq)

fire.seq.lt <- rep(c(rep(0, 20), rep(1, 20), rep(2, 20), rep(3, 20), rep(4, 20)), 2)
temp.seq.lt <- rep(0, 200)
for.seq.lt <- c(rep(-0.8360172, 100), rep(1.1575623, 100))
b50.seq.lt <- rep(c(-1.016001, 0.952501), 100)

######  Models: Effects of long-term fire history on Pr dynamics ########

##### A) Does fire history influence Pr occurrence? ######

preds <- Plot[c("burned1950.c", "tmax.c", "forest.c")]
data <- list(Pr = Plot$Pr2006, N = nrow(Plot), X = preds, I = ncol(preds),
             Nsims=nsims, firesim=b50.seq, forsim=for.seq)
occurrence.b50.06 <- stan(file = "occurrence.stan",
                          data = data, 
                          chains=4 , iter=2000 , warmup=1000, thin=1,
                          control=list(adapt_delta = 0.99, stepsize = 0.01))
#launch_shinystan(occurrence.b50.06)
plot(occurrence.b50.06, ci_level = 0.5, outer_level=0.90, pars=c("a", "betas"))


preds <- Plot[c("totalfire.c", "tmax.c", "forest.c")]
data <- list(Pr = Plot$Pr2006, N = nrow(Plot), X = preds, I = ncol(preds),
             Nsims=nsims, firesim=fire.seq, forsim=for.seq)
occurrence.fires.06 <- stan(file = "occurrence.stan",
                          data = data, 
                          chains=4 , iter=2000 , warmup=1000, thin=1,
                          control=list(adapt_delta = 0.99, stepsize = 0.01))

plot(occurrence.fires.06, ci_level = 0.5, outer_level=0.90, pars=c("a", "betas"))
#launch_shinystan(occurrence.fires.06)

## Posterior predictive checks & Model Fit
y <- Plot$Pr2006
ypred <- rstan::extract(occurrence.fires.06)
ypred <- ypred$ypred

prop_zero <- function(x) mean(x == 1)
ppc_stat(y, ypred, stat = "prop_zero") + ylab("Frequency") + xlab("Proportion of successes/trials") + ggtitle("Model = Occurrence of P. ramorum, with number of fires as variable")

ypred <- rstan::extract(occurrence.b50.06)
ypred <- ypred$ypred
ppc_stat(y, ypred, stat = "prop_zero") + ylab("Frequency") + xlab("Proportion of successes/trials") + ggtitle("Model = Occurrence of P. ramorum, \n with fire history since 1950 as variable")


#Area under the receiver-operating curve
preds <- summary(occurrence.fires.06, pars='p', probs = c(0.25,0.5,0.75))
preds <- as.data.frame(preds$summary)
roc(Plot$Pr2006, plogis(preds$`50%`)) ## 0.70

preds <- summary(occurrence.b50.06, pars='p', probs = c(0.25,0.5,0.75))
preds <- as.data.frame(preds$summary)
roc(Plot$Pr2006, plogis(preds$`50%`)) ## 0.71


##### B) Does fire history influence intensity of Pr infestation? ######
##### Subsetted to just plots where Pr is present, to focus on effects of local disease dynamics, not effects of differences in time since invasion:

Xsim <- data.frame(b50.seq.lt, temp.seq.lt, for.seq)

## effect of burned/unburned since 1950:
preds <- LTPlot[c("burned1950.c", "tmax.c", "forest.c")]
data <- list(infected = LTPlot$totinf13, 
             totalhosts = LTPlot$tothost13, 
             N = nrow(LTPlot), 
             X = preds, 
             I = ncol(preds),
             Nsims=nsims, Xsim=Xsim)

infections.13.b50 <- stan(file = "infectionintensity.stan",
                      data = data, 
                      chains=4 , iter=2000 , warmup=1000, thin=1,
                      control=list(adapt_delta = 0.99, stepsize = 0.01))
plot(infections.13.b50, ci_level = 0.5, outer_level=0.90, pars=c("betas"))


## Posterior predictive checks & Model Fit
y <- LTPlot$totinf13
ypred <- rstan::extract(infections.13.b50)
ypred <- ypred$ypred
ppc_dens_overlay(y, ypred[1:50,]) + ylab("Frequency") + xlab("# of infections observed") + ggtitle("Model = Infestation intensity, \nwith fire history since 1950 as variable")
ypred <- apply(ypred, 2, median)
mae(y/LTPlot$tothost13, ypred/LTPlot$tothost13)

## effect of # of fires:
Xsim <- data.frame(fire.seq.lt, temp.seq.lt, for.seq)
preds <- LTPlot[c("totalfire.c", "tmax.c", "forest.c")]
data <- list(infected = LTPlot$totinf13, 
             totalhosts = LTPlot$tothost13, 
             N = nrow(LTPlot), 
             X = preds, 
             I = ncol(preds), 
             Nsims=nsims, Xsim=Xsim)

infections.13.fires <- stan(file = "infectionintensity.stan",
                          data = data, 
                          chains=4 , iter=2000 , warmup=1000, thin=1,
                          control=list(adapt_delta = 0.99, stepsize = 0.01))

plot(infections.13.fires, ci_level = 0.5, outer_level=0.95, pars=c("a", "betas"))


## Posterior predictive checks & Model Fit
y <- LTPlot$totinf13
ypred <- rstan::extract(infections.13.fires)
ypred <- ypred$ypred
ppc_dens_overlay(y, ypred[1:50,]) + ylab("Frequency") + xlab("# of infections observed") + ggtitle("Model = Infestation intensity, \nwith number of fires since 1950 as variable")

ypred <- apply(ypred, 2, median)
mae(y/LTPlot$tothost13, ypred/LTPlot$tothost13)



# C) Does fire history influence patterns of host mortality in infested areas? ######

LTPlot <- subset(LTPlot, P_ramorum==TRUE & Burned_2008==FALSE) ## Subset to just longterm plots that didn't burn in 2008 (so we can look at fire history, not recent fire effects on tree mortality)

# Rescale variables:
LTPlot$tmax.c <- scale(LTPlot$tmax)
LTPlot$plotindex <- group_indices(LTPlot, BigSurPlot)
LTPlot$ForestAllianceType <- group_indices(LTPlot, ForestAllianceType)
LTPlot$totalfire.c <- scale(LTPlot$totalfire)
LTPlot$burned1950.c <- scale(LTPlot$burned1950.08)
LTPlot$forest.c <- scale(LTPlot$ForestAllianceType)
LTPlot$nonzero.hmort13 <- ifelse(LTPlot$hostmortality.pr>0, 1, 0) ## SOD-related Mortality occurred

# Effect of burned since 1950:
mortality13.b50.z <- stan_glm(nonzero.hmort13 ~ burned1950.c + tmax.c + (forest.c), 
                            data = subset(LTPlot), family = binomial(link = logit))

mortality13.b50.nz <-  stan_glm(hostmortality.pr ~ burned1950.c + tmax.c + (forest.c), 
                              data = subset(LTPlot, nonzero.hmort13==1),  
                              family = Gamma(link= log), 
                              prior = normal(0,1), prior_intercept = normal(0, 10),  
                              chains = 3, cores = 2000, seed = 500, adapt_delta=0.99)

plot(mortality13.b50.z)
plot(mortality13.b50.nz)


## Posterior predictive checks & Model Fit
yrep <- posterior_predict(mortality13.b50.nz, draws = 500)
y <- LTPlot$hostmortality.pr[LTPlot$nonzero.hmort13==1]
ppc_dens_overlay(y, yrep[1:30, ])+ ylab("Density") + xlab("Basal area of host mortality observed") + ggtitle("Model = Host mortality, \nwith fire history since 1950 as variable") 
yrep <- apply(yrep, 2, median)
mae(y, yrep)



## Effect of number of fires since 1950
mortality13.fires.z <- stan_glm(nonzero.hmort13 ~ totalfire.c + tmax.c + (forest.c), 
                              data = subset(LTPlot), 
                              family = binomial(link = logit))

mortality13.fires.nz <-  stan_glm(hostmortality.pr ~ totalfire +tmax.c + (forest.c),
                                data = subset(LTPlot, nonzero.hmort13==1),  
                                family = Gamma(link= log), 
                                prior = normal(0,1), prior_intercept = normal(0, 10),  
                                chains = 3, cores = 2000, seed = 500, adapt_delta=0.99)

plot(mortality13.fires.z)
plot(mortality13.fires.nz)


## Posterior predictive checks & Model Fit
yrep <- posterior_predict(mortality13.fires.nz, draws = 500)
y <- LTPlot$hostmortality.pr[LTPlot$nonzero.hmort13==1]
ppc_dens_overlay(y, yrep[1:30, ])+ ylab("Density") + xlab("Basal area of host mortality observed") + ggtitle("Model = Host mortality, \nwith number of fires since 1950 as variable") 
yrep <- apply(yrep, 2, median)
mae(y, yrep)






