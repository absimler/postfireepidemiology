## Code for "Wildfire alters the disturbance impacts of an emerging forest disease via changes to host occurrence and demographic structure" ###

## Part III. Impacts of recent fire on SOD dynamics #####

library(rstan)
library(rstanarm)
library(bayesplot)
library(loo)
library(tidyverse)


# Load data:
Trees = read.csv("datasets/hosttrees.csv", header=T, na.strings=c("", " "))
Stems = read.csv("datasets/hoststems.csv", header=T, na.strings=c("", " "))
LTPlot <- read.csv("datasets/longtermPlots_201014.csv", header=T, na.strings=c("", " "))


## Add new variables:
LTPlot <- subset(LTPlot, BigSurPlot!="BS187") # Exclude 1 plot within fire perimeter that did not burn
LTPlot$Pr2010 <- ifelse(LTPlot$totinf10>0, 1, 0)
LTPlot$Pr2013 <- ifelse(LTPlot$totinf13>0, 1, 0)
LTPlot$RShosts <- LTPlot$UMCAnew + LTPlot$LIDEnew ## Total resprouting hosts



## Subset unburned data: ####
LTPlot.UB <- subset(LTPlot, Burned_2008==FALSE)
LTPlot.UB$tmax.c <- scale(LTPlot.UB$tmax)
LTPlot.UB$Pr2010.c <- scale(LTPlot.UB$Pr2010)
LTPlot.UB$UMCA13.c <- scale(LTPlot.UB$totinfumca13)
LTPlot.UB$plotindex <- group_indices(LTPlot.UB, BigSurPlot)
LTPlot.UB$ForestAllianceType.c <- scale(group_indices(LTPlot.UB, ForestAllianceType))

## Subset to only stems and trees within those burned plots
Stems <-  Stems %>% drop_na(BA.2010)
LTstems.UB =subset(Stems, Plot %in% unique(LTPlot.UB$BigSurPlot)) 
LTtrees.UB =subset(Trees, Plot %in% unique(LTPlot.UB$BigSurPlot)) 
LTstems.UB$Species <- as.factor(LTstems.UB$Species)
LTstems.UB$BA.2010 <- as.numeric(as.character(LTstems.UB$BA.2010))
LTstems.UB$DBH2010c <- scale(as.numeric(LTstems.UB$DBH.2010))
LTtrees.UB$Prpos1011c <- scale(as.numeric(LTtrees.UB$Prpos1011))

## subset to only include most common susceptible host species:
LTstems.UB <- subset(LTstems.UB, Species=="QUAG" | Species=="LIDE")
LTstems.UB <- subset(LTstems.UB, Status.2010=="L" & (Mort13==1 | Mort13==0) & BA.2010>0)
LTtrees.UB <- subset(LTtrees.UB, Species=="QUAG" | Species=="LIDE")
LTPlot.UB =subset(LTPlot.UB, BigSurPlot %in% unique(LTtrees.UB$Plot)) 
LTtrees.UB =subset(LTtrees.UB, TreeID %in% unique(LTstems.UB$TreeID)) 

LTstems.UB$treeindex <- group_indices(LTstems.UB, TreeID)
LTstems.UB$plotindex <- group_indices(LTstems.UB, Plot)
LTtrees.UB$treeindex <- group_indices(LTtrees.UB, TreeID)
LTtrees.UB$plotindex <- group_indices(LTtrees.UB, Plot)
LTtrees.UB$speciesindex <- group_indices(LTtrees.UB, Species)
LTstems.UB$Mort13 <- as.numeric(as.character(LTstems.UB$Mort13))
LTPlot.UB$burned1950.c <- scale(LTPlot.UB$burned1950)

## Subset burned data ####
LTPlot.B <- subset(LTPlot, Burned_2008==TRUE)
LTPlot.B$tmax.c <- scale(LTPlot.B$tmax)
LTPlot.B$UMCA13.c <- scale(LTPlot.B$totinfumca13)
LTPlot.B$plotindex <- group_indices(LTPlot.B, BigSurPlot)
LTPlot.B$ForestAllianceType.c <- scale(group_indices(LTPlot.B, ForestAllianceType))
LTPlot.B$LiveBA10.c <- scale(LTPlot.B$LiveBA10)
LTPlot.B$LIDEsurvived.c <- scale(LTPlot.B$LIDEsurvived)
LTPlot.B$LIDEnew.c <- scale(LTPlot.B$LIDEnew)
LTPlot.B$UMCAsurvived.c <- scale(LTPlot.B$UMCAsurvived)
LTPlot.B$UMCAnew.c <- scale(LTPlot.B$UMCAnew)
LTPlot.B$Pr2010.c <- scale(LTPlot.B$Pr2010)
LTPlot.B$burned1950.c <- scale(LTPlot.B$burned1950)
LTPlot.B$RShosts.c <- scale(LTPlot.B$RShosts)

## Subset to only trees and stems present in longterm plots
LTstems.B =subset(Stems, Plot %in% unique(LTPlot.B$BigSurPlot)) 
LTtrees.B =subset(Trees, Plot %in% unique(LTPlot.B$BigSurPlot)) 
LTstems.B$Species <- as.factor(LTstems.B$Species)
LTstems.B$BA.2010 <- as.numeric(as.character(LTstems.B$BA.2010))
LTstems.B$DBH2010c <- scale(as.numeric(LTstems.B$DBH.2010))
LTtrees.B$Prpos1011c <- scale(as.numeric(LTtrees.B$Prpos1011))

## For the purposes of this analysis, we only care about susceptible hosts (LIDE, QUER spp), so sBset to only include those, and only the most common susceptible host species:
LTstems.B <- subset(LTstems.B, Species=="QUAG" | Species=="LIDE")
LTstems.B <- subset(LTstems.B, Status.2010=="L" & (Mort13==1 | Mort13==0) & BA.2010>0)
LTtrees.B <- subset(LTtrees.B, Species=="QUAG" |  Species=="LIDE")
LTPlot.B =subset(LTPlot.B, BigSurPlot %in% unique(LTtrees.B$Plot)) 
LTtrees.B =subset(LTtrees.B, TreeID %in% unique(LTstems.B$TreeID)) 

LTstems.B$treeindex <- group_indices(LTstems.B, TreeID)
LTstems.B$plotindex <- group_indices(LTstems.B, Plot)
LTPlot.B$plotindex <- group_indices(LTPlot.B, BigSurPlot)
LTtrees.B$plotindex <- group_indices(LTtrees.B, Plot)
LTtrees.B$speciesindex <- group_indices(LTtrees.B, Species)
LTstems.B$Mort13 <- as.numeric(as.character(LTstems.B$Mort13))


#### A) How does fire influence the effects of P. ramorum on stem mortality? ################

## unburned model:
plotpred <- LTPlot.UB[c("ForestAllianceType.c", "tmax.c",  
                        "burned1950.c", "UMCA13.c")]
treepred <- LTtrees.UB[c("speciesindex", "Prpos1011c")]
stempred <- LTstems.UB[c("DBH2010c")]

data <- list( S = nrow(LTstems.UB), 
              TR = nrow(LTtrees.UB), 
              P = nrow(LTPlot.UB), 
              sx = stempred, 
              px = plotpred, 
              tx = treepred,
              TreeID = LTstems.UB$treeindex, 
              PlotID=LTstems.UB$plotindex, 
              TreeinPlot = LTtrees.UB$plotindex, 
              J = ncol(stempred), 
              K = ncol(treepred),
              L = ncol(plotpred),
              Dead=LTstems.UB$Mort13)
              
stemmort.unburned <- stan(file = "stemmortality.stan",
                          data = data, 
                      chains=2 , iter=2000 , warmup=1000,
                      control=list(adapt_delta = 0.999, max_treedepth=20))

plot(stemmort.unburned, ci_level = 0.5, outer_level=0.90, 
     pars=c("mu_alpha","sbeta","tbeta", "pbeta"))


### Burned model:
plotpred <- LTPlot.B[c("ForestAllianceType.c", "tmax.c","UMCA13.c")]
treepred <- LTtrees.B[c("speciesindex", "Prpos1011c")]
stempred <- LTstems.B[c("DBH2010c")]


data <- list( S = nrow(LTstems.B), 
              TR = nrow(LTtrees.B), 
              P = nrow(LTPlot.B), 
              sx = stempred, 
              px = plotpred, 
              tx = treepred,
              TreeID = LTstems.B$treeindex, 
              PlotID=LTstems.B$plotindex, 
              TreeinPlot = LTtrees.B$plotindex, 
              J = ncol(stempred), 
              K = ncol(treepred),
              L = ncol(plotpred),
              Dead=LTstems.B$Mort13)

stemmort.burned <- stan(file = "stemmortality.stan",
                          data = data, 
                          chains=2 , iter=2000 , warmup=1000, thin=1,
                          control=list(adapt_delta = 0.99, max_treedepth=20))

plot(stemmort.burned, ci_level = 0.5, outer_level=0.90, 
     pars=c("mu_alpha", "sbeta", "tbeta", "pbeta"))


## Null model burned

data <- list( S = nrow(LTstems.B), 
              TR = nrow(LTtrees.B), 
              P = nrow(LTPlot.B), 
              TreeID = LTstems.B$treeindex, 
              PlotID=LTstems.B$plotindex, 
              TreeinPlot = LTtrees.B$plotindex, 
              Dead=LTstems.B$Mort13)

stemmort.burned.null <- stan(file = "stemmortalitynull.stan",
                        data = data, 
                        chains=2 , iter=2000 , warmup=500, thin=1,
                        control=list(adapt_delta = 0.99, max_treedepth=20))

## Null model unburned

data <- list( S = nrow(LTstems.UB), 
              TR = nrow(LTtrees.UB), 
              P = nrow(LTPlot.UB), 
              TreeID = LTstems.UB$treeindex, 
              PlotID=LTstems.UB$plotindex, 
              TreeinPlot = LTtrees.UB$plotindex, 
              Dead=LTstems.UB$Mort13)

stemmort.unburned.null <- stan(file = "stemmortalitynull.stan",
                             data = data, 
                             chains=2 , iter=2000 , warmup=500, thin=1,
                             control=list(adapt_delta = 0.99, max_treedepth=20))


####Comparison of looic between null and full models:

loo.mortburned <- loo(extract_log_lik(stemmort.burned, parameter_name = "log_lik"))
loo.mortunburned <- loo(extract_log_lik(stemmort.unburned, parameter_name = "log_lik"))
loo.mortburned.null <- loo(extract_log_lik(stemmort.burned.null, parameter_name = "log_lik"))
loo.mortunburned.null <- loo(extract_log_lik(stemmort.unburned.null, parameter_name = "log_lik"))

compare(loo.mortburned.null, loo.mortburned)
compare(loo.mortunburned.null, loo.mortunburned)


## B) Which factors determine intensity of post-fire infestation? ####################

## Subset data to just burned plots where P. ramorum has been confirmed & where hosts are present post-fire (don't want to model invasion, but rather how disease proceeds in regions where it is known to occur)
LTPlot.B2 <- subset(LTPlot,  P_ramorum==TRUE & Burned_2008==TRUE & hostpresent13>0)

LTPlot.B2$tmax.c <- scale(LTPlot.B2$tmax)
LTPlot.B2$UMCABA13.c <- scale(LTPlot.B2$UMCA.BA.LIVE.13)
LTPlot.B2$plotindex <- group_indices(LTPlot.B2, BigSurPlot)
LTPlot.B2$ForestAllianceType.c <- scale(group_indices(LTPlot.B2, ForestAllianceType))
LTPlot.B2$LiveBA10.c <- scale(LTPlot.B2$LiveBA10)
LTPlot.B2$LIDEsurvived.c <- scale(LTPlot.B2$LIDEsurvived)
LTPlot.B2$LIDEnew.c <- scale(LTPlot.B2$LIDEnew)
LTPlot.B2$UMCAsurvived.c <- scale(LTPlot.B2$UMCAsurvived)
LTPlot.B2$RShosts.c <- scale(LTPlot.B2$RShosts)
LTPlot.B2$Pr2010.c <- scale(LTPlot.B2$Pr2010)
LTPlot.B2$burned1950.c <- scale(LTPlot.B2$burned1950)


## Simulation data for plotting marginal effects
UMCA.seq <- rep(seq(from=min(LTPlot.B2$UMCAsurvived.c), 
                    to=max(LTPlot.B2$UMCAsurvived.c), length.out = 100), 2)
for.seq <- c(rep(-0.884, 100), rep(1.106, 100))
tmax.seq <- rep(mean(LTPlot.B2$tmax.c))
BA.seq <- rep(mean(LTPlot.B2$LiveBA10.c))
LIDE.seq <- rep(mean(LTPlot.B2$LIDEsurvived.c))
LIDERS.seq <- rep(0, 200)
UMCARS.seq <- rep(0, 200)
RS.seq <- rep(mean(LTPlot.B2$RShosts.c))
nsims <- length(for.seq)



preds <- LTPlot.B2[c( "ForestAllianceType.c", "tmax.c",  "LiveBA10.c",
                      "LIDEsurvived.c", "UMCAsurvived.c", "RShosts.c")]

Xseq <- cbind(for.seq, tmax.seq, BA.seq, LIDE.seq, UMCA.seq, RS.seq)


data <- list(infected = LTPlot.B2$totinf13, 
             totalhosts = LTPlot.B2$tothost13,
             N = nrow(LTPlot.B2), 
             X = preds, 
             I = ncol(preds), 
             Nsims=nsims, Xsim = Xseq)

infections.13 <- stan(file = "infectionintensity.stan",
                          data = data, 
                          chains=3 , iter=2000 , warmup=500, thin=1,
                          control=list(adapt_delta = 0.99, stepsize = 0.01))

plot(infections.13, ci_level = 0.5, outer_level=0.95, pars=c("a", "betas"))

y <- LTPlot.B2$totinf13
ypred <- rstan::extract(infections.13)
ypred <- ypred$ypred
ppc_dens_overlay(y, ypred[1:50,]) + ylab("Density") + xlab("# of infections observed") + ggtitle("Model = Post-fire infection of all hosts")
ypred <- apply(ypred, 2, median)
mae(y/LTPlot.B2$tothost13, ypred/LTPlot.B2$tothost13)


# Null model:
data <- list(infected = LTPlot.B2$totinf13, totalhosts = LTPlot.B2$tothost13,
             N = nrow(LTPlot.B2),
             Nsims=nsims)
infections.null <- stan(file = "infectionintensitynull.stan",
                      data = data, 
                      chains=3 , iter=2000 , warmup=1000, thin=1,
                      control=list(adapt_delta = 0.99, stepsize = 0.01))

loglik <- extract_log_lik(infections.null, parameter_name = "log_lik")
loo.infectionsnull <- loo(loglik, threshold=0.7)
loglik <- extract_log_lik(infections.13, parameter_name = "log_lik")
loo.infections13 <- loo(loglik, threshold=0.7)

compare(loo.infectionsnull, loo.infections13) 



## C) Which factors determine post-fire re-infection? #######

## Simulation data for plotting marginal effects:
UMCA.seq <- rep(seq(from=min(LTPlot.B2$UMCAsurvived.c), 
                    to=max(LTPlot.B2$UMCAsurvived.c), length.out = 100), 2)
for.seq <- c(rep(-0.884, 100), rep(1.106, 100))
tmax.seq <- rep(mean(LTPlot.B2$tmax.c))
BA.seq <- rep(mean(LTPlot.B2$LiveBA10.c))
LIDE.seq <- rep(mean(LTPlot.B2$LIDEsurvived.c))
RS.seq <- rep(mean(LTPlot.B2$RShosts.c))

nsims <- length(for.seq)
preds <- LTPlot.B2[c( "ForestAllianceType.c", "tmax.c",  "LiveBA10.c",
                     "LIDEsurvived.c", "UMCAsurvived.c", "RShosts.c")]

Xseq <- cbind(for.seq, tmax.seq, BA.seq, LIDE.seq, UMCA.seq, RS.seq)

data <- list(infected = LTPlot.B2$topinf13, 
             totalhosts = LTPlot.B2$topkillhost13,
             N = nrow(LTPlot.B2), 
             X = preds, 
             I = ncol(preds), 
             Nsims=nsims, Xsim = Xseq)

topkillinf.13 <- stan(file = "infectionintensity.stan",
                      data = data, 
                      chains=3 , iter=2000 , warmup=1000, thin=1,
                      control=list(adapt_delta = 0.99, stepsize = 0.01))

plot(topkillinf.13, ci_level = 0.5, outer_level=0.90, pars=c("a", "betas"))

y <- LTPlot.B2$topinf13
ypred <- rstan::extract(topkillinf.13)
ypred <- ypred$ypred
ppc_dens_overlay(y, ypred[1:50,]) + ylab("Density") + xlab("# of infections observed") + ggtitle("Model = Post-fire infection of topkilled hosts")
ypred <- apply(ypred, 2, median)
mae(y/LTPlot.B2$topkillhost13, ypred/LTPlot.B2$topkillhost13)



# Null model:
data <- list(infected = LTPlot.B2$topinf13, totalhosts = LTPlot.B2$topkillhost13,
             N = nrow(LTPlot.B2),
             Nsims=nsims, Xsim = Xseq)
topkillinf.null <- stan(file = "infectionintensitynull.stan",
                        data = data, 
                        chains=3 , iter=2000 , warmup=1000, thin=1,
                        control=list(adapt_delta = 0.99, stepsize = 0.01))

loglik <- extract_log_lik(topkillinf.null, parameter_name = "log_lik")
loo.infectionsnullTK <- loo(loglik, threshold=0.7)
loglik <- extract_log_lik(topkillinf.13, parameter_name = "log_lik")
loo.infections13TK <- loo(loglik, threshold=0.7)

compare(loo.infectionsnullTK, loo.infections13TK)
