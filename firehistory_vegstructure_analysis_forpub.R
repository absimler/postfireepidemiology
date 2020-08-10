## Code for "Wildfire alters the disturbance impacts of an emerging forest disease via changes to host occurrence and demographic structure" ###

### Part I. Do key aspects of forest structure & composition vary with fire history? #####

library(rstanarm)
library(rstan)
library(tidyverse)
library(pROC)
library(loo)

## Load set of all 280 plots in Big Sur Network:
Plot = read.csv("datasets/all280Plots_20062007.csv", header=T, na.strings=c("", " "))

## Create categorical variables for presence/absence of hosts:
Plot$UMCAnonzero <- ifelse(Plot$UMCA.BA.LIVE.0607>0, 1, 0)
Plot$LIDEnonzero <- ifelse(Plot$LIDE.BA.LIVE.0607>0, 1, 0)
Plot$HOSTnonzero <- ifelse(Plot$Hosts06>0, 1, 0)

## Scale and center variables:
Plot$tmax <- scale(Plot$tmax)
Plot$totalfire.nb <- scale(Plot$totalfire.nb)

# MAE function
mae <- function(y, ypred) {
  mae = (mean(abs(y - ypred)))
  return(mae)
}

#r2 function
r2<-function(y_hat,y) {
    RSS<-sum((((y_hat))-(y))^2)
    TSS<-sum(((y)-(mean(y)))^2)
    return ( 1 -RSS/TSS)}



## A) Effect of fire history on plot-level basal area ##########################

## simulation data for plotting marginal effects:
fire.seq <- rep(c(rep(0, 20), rep(1, 20), rep(2, 20), rep(3, 20), rep(4, 20)), 2)
b50.seq <- rep(c(0,1), 100)
for.seq <- c(rep(1, 100), rep(2, 100))
nsims <- length(fire.seq)
temp.seq <- rep(mean(Plot$tmax), 200)

# A1) Variable for total number of burns since 1950
preds <- Plot[c( "ForestAllianceType", "tmax",  "totalfire.nb")]
Xsim = cbind(for.seq, temp.seq, fire.seq)

data <- list(y= Plot$LiveBA0607, X = preds, N = nrow(Plot),
             Xsim = Xsim, I = ncol(preds),
             Nsim= nsims)

totalBA.tot <- stan(file = "half_normal2.stan",
                    data = data, 
                    chains=4 , iter=2000 , warmup=1000, thin=1,
                    control=list(adapt_delta = 0.99, stepsize = 0.01))

plot(totalBA.tot, outer_level=0.9) # no effect

fit <- rstan::extract(totalBA.tot)
yhat <- apply(fit$mu, 2, median) 

mae(Plot$LiveBA0607, yhat)
r2(yhat, Plot$LiveBA0607)

# A2) Variable for whether plot burned since 1950
preds <- Plot[c( "ForestAllianceType", "tmax",  "burned1950")]
Xsim = cbind(for.seq, temp.seq, b50.seq)

data <- list(y= Plot$LiveBA0607, X = preds, N = nrow(Plot),
             Xsim = Xsim, I = ncol(preds),
             Nsim= nsims)
totalBA.burned <- stan(file = "half_normal2.stan",
                    data = data, 
                    chains=4 , iter=2000 , warmup=1000, thin=1,
                    control=list(adapt_delta = 0.99, stepsize = 0.01))

plot(totalBA.burned, pars="betas", outer_level=0.9) 

fit <- rstan::extract(totalBA.burned)
yhat <- apply(fit$mu, 2, median) 

mae(Plot$LiveBA0607, yhat)
r2(yhat, Plot$LiveBA0607)


## B) Effect of fire history on plot-level basal area of SOD hosts ########################

# 1B) Variable for Total number of burns since 1950
preds <- Plot[c( "ForestAllianceType", "tmax",  "totalfire.nb")]
Xsim = cbind(for.seq, temp.seq, fire.seq)

data <- list(y= Plot$Hosts06, X = preds, N = nrow(Plot),
             Xsim = Xsim, I = ncol(preds),
             Nsim= nsims)

hostBA.tot <- stan(file = "half_normal2.stan",
                    data = data, 
                    chains=4 , iter=2000 , warmup=1000, thin=1,
                    control=list(adapt_delta = 0.99, stepsize = 0.01))

plot(hostBA.tot, pars="betas", outer_level=0.9) 

fit <- rstan::extract(hostBA.tot)
yhat <- apply(fit$mu, 2, median) 

mae(Plot$Hosts06, yhat)
r2(yhat, Plot$Hosts06)


# B2) Variable for whether plot burned since 1950
preds <- Plot[c("ForestAllianceType", "tmax",  "burned1950")]
Xsim = cbind(for.seq, temp.seq, b50.seq)

data <- list(y= Plot$Hosts06, X = preds, N = nrow(Plot),
             Xsim = Xsim, I = ncol(preds),
             Nsim= nsims)
hostBA.burned <- stan(file = "half_normal2.stan",
                       data = data, 
                       chains=4 , iter=2000 , warmup=1000, thin=1,
                       control=list(adapt_delta = 0.99, stepsize = 0.01))

plot(hostBA.burned, pars="betas", outer_level=0.9)

fit <- rstan::extract(hostBA.burned)
yhat <- apply(fit$mu, 2, median) 

mae(Plot$Hosts06, yhat)
r2(yhat, Plot$Hosts06)



## C) Effect of fire history on BAY basal area in plots: Hurdle approach #########

## C1) Plot burned since 1950:
Bay.burned1950.zero <- stan_glm(UMCAnonzero ~ burned1950 + 
                                  (ForestAllianceType) + tmax, 
                                  data = Plot,  
                                family = binomial(link = logit), 
                                  prior = normal(0,1), 
                                prior_intercept = normal(0, 10),  
                                  chains = 3, cores = 2000, 
                                seed = 500, adapt_delta=0.99)

plot(Bay.burned1950.zero, pars="burned1950")

mu <- posterior_linpred(Bay.burned1950.zero, transform = TRUE)
mu <- apply(mu, 2, median) 
roc(Plot$UMCAnonzero, mu)

Bay.burned1950.nonzero <- stan_glm(UMCA.BA.LIVE.0607 ~ burned1950 + 
                                     (ForestAllianceType) + tmax, 
                                     data = subset(Plot, UMCAnonzero==1),  
                                   family = Gamma(link = log), 
                                     prior = normal(0,1), 
                                   prior_intercept = normal(0, 10),  
                                     chains = 3, cores = 2000, 
                                   seed = 500, adapt_delta=0.99)

plot(Bay.burned1950.nonzero, pars="burned1950")

mu <- posterior_predict(Bay.burned1950.nonzero, transform = TRUE)
mu <- apply(mu, 2, mean) 
mae(Plot$UMCA.BA.LIVE.0607[Plot$UMCAnonzero==1], mu)


# C2) Number of fires since 1950
Bay.totalfire.zero <- stan_glm(UMCAnonzero ~ totalfire.nb + tmax +
                                 (ForestAllianceType), 
                                 data = Plot,  family = binomial(link = logit), 
                                 prior = normal(0,1), prior_intercept = normal(0, 10),  
                                 chains = 3, cores = 2000, seed = 500, adapt_delta=0.99)

plot(Bay.totalfire.zero)

mu <- posterior_linpred(Bay.totalfire.zero, transform = TRUE)
mu <- apply(mu, 2, median) 
roc(Plot$UMCAnonzero, mu)


Bay.totalfire.nonzero <- stan_glm(UMCA.BA.LIVE.0607 ~ tmax +
                                    totalfire.nb + (ForestAllianceType), 
                                    data = subset(Plot, UMCAnonzero==1),  
                                  family = Gamma(link = log),
                                    prior = normal(0,1), 
                                  prior_intercept = normal(0, 10),  
                                    chains = 3, cores = 2000, 
                                  seed = 500, adapt_delta=0.99)

mu <- posterior_predict(Bay.totalfire.nonzero, transform = TRUE)
mu <- apply(mu, 2, mean) 
mae(Plot$UMCA.BA.LIVE.0607[Plot$UMCAnonzero==1], mu)

plot(Bay.totalfire.nonzero)

yrep <- posterior_predict(Bay.burned1950.nonzero, draws=55)
y <- Plot$UMCA.BA.LIVE.10[Plot$UMCAnonzero==1]



## D) Effect of fire history on Tanoak basal area in plots: #########

## D1) Plot burned since 1950:
Tanoak.burned1950.zero <- stan_glm(LIDEnonzero ~ burned1950 + tmax +
                                     (ForestAllianceType), 
                                  data = Plot,  family = binomial(link = logit), 
                                  prior = normal(0,1), 
                                  prior_intercept = normal(0, 10),  
                                  chains = 3, cores = 2000, 
                                  seed = 500, adapt_delta=0.99)

plot(Tanoak.burned1950.zero)

mu <- posterior_linpred(Tanoak.burned1950.zero, transform = TRUE)
mu <- apply(mu, 2, median) 
roc(Plot$LIDEnonzero, mu)

Tanoak.burned1950.nonzero <- stan_glm(LIDE.BA.LIVE.0607 ~ burned1950 + tmax +
                                        (ForestAllianceType), 
                                     data = subset(Plot, LIDEnonzero==1),  
                                     family = Gamma(link = log), 
                                     prior = normal(0,1), 
                                     prior_intercept = normal(0, 10),  
                                     chains = 3, cores = 2000, 
                                     seed = 500, adapt_delta=0.99)

plot(Tanoak.burned1950.nonzero)

mu <- posterior_predict(Tanoak.burned1950.nonzero, transform = TRUE)
mu <- apply(mu, 2, mean) 
mae(Plot$LIDE.BA.LIVE.0607[Plot$LIDEnonzero==1], mu)


# D2) number of fires since 1950
Tanoak.totalfire.zero <- stan_glm(LIDEnonzero ~ totalfire.nb + tmax
                                  + (ForestAllianceType), 
                                 data = Plot,  family = binomial(link = logit), 
                                 prior = normal(0,1), prior_intercept = normal(0, 10),  
                                 chains = 3, cores = 2000, seed = 500, adapt_delta=0.99)

plot(Tanoak.totalfire.zero)

mu <- posterior_linpred(Tanoak.totalfire.zero, transform = TRUE)
mu <- apply(mu, 2, median) 
roc(Plot$LIDEnonzero, mu)


Tanoak.totalfire.nonzero <- stan_glm(LIDE.BA.LIVE.0607 ~ totalfire.nb + 
                                        tmax +
                                        (ForestAllianceType), 
                                      data = subset(Plot, LIDEnonzero==1),  
                                      family = Gamma(link = log), 
                                      prior = normal(0,1), 
                                      prior_intercept = normal(0, 10),  
                                      chains = 3, cores = 2000, 
                                      seed = 500, adapt_delta=0.99)

plot(Tanoak.totalfire.nonzero)

mu <- posterior_predict(Tanoak.totalfire.nonzero, transform = TRUE)
mu <- apply(mu, 2, mean) 
mae(Plot$LIDE.BA.LIVE.0607[Plot$LIDEnonzero==1], mu)


