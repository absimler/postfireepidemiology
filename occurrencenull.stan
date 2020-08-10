data{
  int<lower=1> N; //number of trees
  int Pr[N]; //the response variable,Pr presence
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
  Pr ~ binomial_logit( 1 , p ); //likelihood
}
generated quantities{
  vector[N] p;
 vector[N] log_lik; //for calculating loo
  for ( i in 1:N ) {
  p[i] = a ;
  log_lik[i] =binomial_logit_lpmf( Pr[i] | 1 , p[i] );
  }
}
  
