data {
  int<lower=1> N; 
  int<lower=1> I;
  real<lower=0> y[N];
}
parameters {
  real a;
  real<lower=0> sigma_y;
}
model {
  a ~ normal(0, 10);
  sigma_y ~ cauchy(0, 5);
  y ~ normal(a, sigma_y);
}
generated quantities{
real log_lik[N];
for (n in 1:N){
 log_lik[n] = normal_lpdf(y[n]| a, sigma_y);
  }
}