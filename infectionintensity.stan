data{
  int<lower=1> N; // # of plots
  int<lower=0> infected[N]; // # of infected
  int<lower=1> totalhosts[N]; // # of total hosts available
  int<lower=0> I; // # predictors of interest
  matrix[N, I] X; // predictors of interest
  int<lower=1> Nsims;
  matrix[Nsims, I] Xsim; // predictors of interest
}
parameters{
  real a;
  vector[I] betas;
}
model{
  vector[N] p;
  betas ~ normal( 0 , 1);
  a ~ normal( 0 , 5);
  for ( i in 1:N ) {
    p[i] = a + X[i]*betas;
  }
  infected ~ binomial_logit( totalhosts , p );
}
generated quantities{
vector[N] p;
vector[Nsims] p_sims; //simulated probability for plots
vector[N] log_lik; //for calculating loo
  for ( i in 1:N ) {
  p[i] = a + X[i]*betas;
  log_lik[i] =binomial_logit_lpmf( infected[i] | totalhosts[i] , p[i] );
  }
   for (q in 1:Nsims){
    p_sims[q] = inv_logit(a + Xsim[q]*betas);
  }
}
  
