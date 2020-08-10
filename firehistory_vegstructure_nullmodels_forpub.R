## Code for "Wildfire alters the disturbance impacts of an emerging forest disease via changes to host occurrence and demographic structure" ###

## Null models & Loo scores for Part I:

## Use this file after the Part I/fire history and vegetation structure analysis datasets have been loaded and models run.


#Null model for basal area

data <- list(y= Plot$LiveBA0607, N = nrow(Plot), I = ncol(preds),
             Nsim= nsims)
totalBA.null <- stan(file = "half_normal_null.stan",
                     data = data, 
                     chains=4 , iter=2000 , warmup=1000, thin=1,
                     control=list(adapt_delta = 0.99, stepsize = 0.01))

#Null model for host basal area
data <- list(y= Plot$Hosts06,  N = nrow(Plot),
             I = ncol(preds),
             Nsim= nsims)
hostBA.null <- stan(file = "half_normal_null.stan",
                    data = data, 
                    chains=4 , iter=2000 , warmup=1000, thin=1,
                    control=list(adapt_delta = 0.99, stepsize = 0.01))




## Null models for tanoak and bay presence/abundance:
tanoak.zero.null <- stan_glm(LIDEnonzero ~ 1, 
                             data = Plot,  family = binomial(link = logit), 
                             prior = normal(0,1), prior_intercept = normal(0, 10),  
                             chains = 4, iter = 2000, seed = 500, adapt_delta=0.99)

bay.zero.null <- stan_glm(UMCAnonzero ~ 1, 
                          data = Plot,  family = binomial(link = logit), 
                          prior = normal(0,1), prior_intercept = normal(0, 10),  
                          chains = 4, iter = 2000, seed = 500, adapt_delta=0.99)

Bay.nonzero.null <- stan_glm(UMCA.BA.LIVE.0607 ~ 1, 
                             data = subset(Plot, UMCAnonzero==1),  
                             family = Gamma(link = log), 
                             prior = normal(0,1), 
                             prior_intercept = normal(0, 10),  
                             chains = 4, iter = 2000, 
                             seed = 500, adapt_delta=0.99)

tanoak.nonzero.null <- stan_glm(LIDE.BA.LIVE.0607 ~ 1, 
                                data = subset(Plot, LIDEnonzero==1),  
                                family = Gamma(link = log), 
                                prior = normal(0,1), 
                                prior_intercept = normal(0, 10),  
                                chains = 3, iter = 2000, 
                                seed = 500, adapt_delta=0.99)


## Extract log likelihoods and calculate LOOIC for each model and associated null.
loglik <- extract_log_lik(totalBA.burned, parameter_name = "log_lik",
                         merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))
loo.totalBA <- loo(loglik, r_eff=r_eff, k_threshold=0.7)

loglik <- extract_log_lik(totalBA.tot, parameter_name = "log_lik",     
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))
loo.totalBA.tot <- loo(loglik, r_eff=r_eff, k_threshold=0.7)

loglik <- extract_log_lik(totalBA.null, parameter_name = "log_lik",     
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))                  
loo.BA.null <- loo(loglik, r_eff=r_eff, k_threshold=0.7)

loglik <- extract_log_lik(hostBA.tot, parameter_name = "log_lik", 
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))   
loo.hosts.tot <- loo(loglik, r_eff=r_eff, k_threshold=0.7)

loglik <- extract_log_lik(hostBA.burned, parameter_name = "log_lik", 
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))   
loo.hosts.b50 <- loo(loglik, r_eff=r_eff, k_threshold=0.7)

loglik <- extract_log_lik(hostBA.null, parameter_name = "log_lik", 
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))   
loo.hosts.null <- loo(loglik,r_eff=r_eff, k_threshold=0.7)

loo.Bayb50.nonzero <- loo(Bay.burned1950.nonzero)
loo.Bayfire.nonzero <- loo(Bay.totalfire.nonzero)
loo.Bay.null <- loo(Bay.nonzero.null)

loo.Bayb50.zero <- loo(Bay.burned1950.zero)
loo.Bayfire.zero <- loo(Bay.totalfire.zero)
loo.Bay.null2 <- loo(bay.zero.null, k_threshold=0.7)

loo.tanoakb50.nonzero <- loo(Tanoak.burned1950.nonzero)
loo.tanoaktotalfire.nonzero <- loo(Tanoak.totalfire.nonzero)
loo.tanoakb50.zero <- loo(Tanoak.burned1950.zero)
loo.tanoaktotalfire.zero <- loo(Tanoak.totalfire.zero)
loo.tanoaknonzero.null <- loo(tanoak.nonzero.null)
loo.tanoakzero.null <- loo(tanoak.zero.null, k_threshold=0.7)


## Comparison of LOOIC:
loo_compare(loo.Bay.null, loo.Bayb50.nonzero)
loo_compare(loo.Bay.null, loo.Bayfire.nonzero)
loo_compare(loo.Bay.null2, loo.Bayb50.zero)
loo_compare(loo.Bay.null2, loo.Bayfire.zero)

loo_compare(loo.tanoaknonzero.null, loo.tanoakb50.nonzero)
loo_compare(loo.tanoaknonzero.null, loo.tanoaktotalfire.nonzero)
loo_compare(loo.tanoakzero.null, loo.tanoakb50.zero)
loo_compare(loo.tanoakzero.null, loo.tanoaktotalfire.zero)

loo_compare(loo.BA.null, loo.totalBA)
loo_compare(loo.hosts.null, loo.hosts.b50)
loo_compare(loo.hosts.null, loo.hosts.tot)


