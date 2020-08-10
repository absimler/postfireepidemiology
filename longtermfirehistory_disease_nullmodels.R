## Code for "Wildfire alters the disturbance impacts of an emerging forest disease via changes to host occurrence and demographic structure" ###

## Null models & Loo scores for Part II/long term fire history and disease analysis.

## Run after running Part II analysis

## Null model for occurrence model
data <- list(Pr = Plot$Pr2006,
             N = nrow(Plot))

occurrence.06.null <- stan(file = "occurrencenull.stan",
                           data = data, 
                           chains=4 , iter=2000 , warmup=1000, thin=1,
                           control=list(adapt_delta = 0.99, stepsize = 0.01))


loglik <- extract_log_lik(occurrence.06.null, parameter_name = "log_lik", 
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))   
loo.null.06 <- loo(loglik, r_eff=r_eff)

loglik <- extract_log_lik(occurrence.b50.06, parameter_name = "log_lik", 
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))  
loo.b50.06occur <- loo(loglik, r_eff=r_eff)

loglik <- extract_log_lik(occurrence.fires.06, parameter_name = "log_lik", 
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))  
loo.fires.06occur <- loo(loglik, r_eff=r_eff)

loo_compare(loo.null.06, loo.fires.06occur)
loo_compare(loo.null.06, loo.b50.06occur)


## Nulls for Infection intensity model:
data <- list(infected = LTPlot$totinf13, 
             totalhosts = LTPlot$tothost13, 
             N = nrow(LTPlot))

infections.13.null <- stan(file = "infectionintensitynull.stan",
                           data = data, 
                           chains=4 , iter=2000 , warmup=1000, thin=1,
                           control=list(adapt_delta = 0.99, stepsize = 0.01))


loglik <- extract_log_lik(infections.13.null, parameter_name = "log_lik", 
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))  
loo.null.13 <- loo(loglik, k_threshold=0.7, r_eff=r_eff)
loglik <- extract_log_lik(infections.13.fires, parameter_name = "log_lik", 
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))  
loo.fires.13 <- loo(loglik, k_threshold=0.7, r_eff=r_eff)
loglik <- extract_log_lik(infections.13.b50, parameter_name = "log_lik", 
                          merge_chains = FALSE)
r_eff <- relative_eff(exp(loglik))  
loo.b50.13 <- loo(loglik, k_threshold=0.7, r_eff=r_eff)

compare(loo.null.13, loo.b50.13)
compare(loo.null.13, loo.fires.13)


## Nulls for Severity of mortality model:
mortality13.null.z <-  stan_glm(nonzero.hmort13 ~ 1, 
                                data = subset(LTPlot), family = binomial(link = logit))
mortality13.null.nz <-  stan_glm(hostmortality.pr ~ 1, 
                                 data = subset(LTPlot, nonzero.hmort13==1),  
                            family = Gamma(link= log), 
                            prior = normal(0,1), prior_intercept = normal(0, 10),  
                                 chains = 4, iter = 2000, seed = 500, adapt_delta=0.99)


mortfire.nz.loo <- loo(mortality13.fires.nz, k_threshold=0.7)
mortfire.z.loo <- loo(mortality13.fires.z, k_threshold=0.7)
mortb50nz.loo <- loo(mortality13.b50.nz, k_threshold=0.7)
mortb50z.loo <- loo(mortality13.b50.z, k_threshold=0.7)
mortz.null <- loo(mortality13.null.z, k_threshold=0.7)
mortnz.null <- loo(mortality13.null.nz, k_threshold=0.7)

compare(mortfire.z.loo, mortz.null)  
compare(mortfire.nz.loo, mortnz.null)

compare(mortb50z.loo, mortz.null) 
compare(mortb50nz.loo, mortnz.null)
