data{
  int<lower=1> N; // # of plots
  int<lower=0> infected[N]; // # of infected
  int<lower=1> totalhosts[N]; // # of total hosts available
}
parameters{
  real a;
}
model{
  vector[N] p;
  a ~ normal( 0 , 5);
  for ( i in 1:N ) {
    p[i] = a ;
  }
  infected ~ binomial_logit( totalhosts , p );
}
generated quantities{
vector[N] p;
vector[N] log_lik; //for calculating loo
  for ( i in 1:N ) {
  p[i] = a ;
  log_lik[i] =binomial_logit_lpmf(infected[i] | totalhosts[i] , p[i] );
  }
}
  
