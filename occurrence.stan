data{
  int<lower=1> N; //number of trees
  int Pr[N]; //the response variable,Pr presence
  int<lower=0> I; // # predictors of interest
  matrix[N, I] X; // predictors of interest
  int<lower=1> Nsims;
  vector[Nsims] firesim;
  vector[Nsims] forsim;
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
  Pr ~ binomial_logit( 1 , p ); //likelihood
}
generated quantities{
  vector[Nsims] p_sims; //simulated probability for plots
  vector[Nsims] y_sims; //simulated outcomes for plotting
  vector[N] p;
 vector[N] log_lik; //for calculating loo
  for ( i in 1:N ) {
  p[i] = a + X[i]*betas;
  log_lik[i] =binomial_logit_lpmf( Pr[i] | 1 , p[i] );
  }
    for (q in 1:Nsims){
    p_sims[q] = inv_logit(a + firesim[q]*betas[1] + forsim[q]*betas[3]);
    y_sims[q] = bernoulli_rng(p_sims[q]);
  }
}
  