data {
  int<lower=1> N; 
  int<lower=1> I;
  int<lower=1> Nsim; //number of sims for plotting
   matrix[N, I] X; // predictors of interest
  real<lower=0> y[N];
  matrix[Nsim, I] Xsim; // simulated data
}
parameters {
  real a;
  vector[I] betas;
  real<lower=0> sigma_y;
}
model {
  vector[N] y_hat;
  y_hat = a + X*betas;
  betas ~ normal(0, 5);
  a ~ normal(0, 10);
  sigma_y ~ cauchy(0, 5);
  y ~ normal(y_hat, sigma_y);
}
generated quantities{
vector[Nsim] y_rep;
real log_lik[N];
vector[N] mu;
y_rep = a + Xsim*betas;
for (n in 1:N){
 mu = a + X*betas;
 log_lik[n] = normal_lpdf(y[n]| mu[n], sigma_y);
  }
}