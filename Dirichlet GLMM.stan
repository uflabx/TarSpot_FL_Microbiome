//============================================//
//        Dirichlet GLMM                      //
//============================================//

data {
  int N;           // Number of observations
  int P;           // Number of fixed effects
  int Q;
  int K;           // Number of response classes
  matrix [N, P] X; // FE design matrix
  matrix [N, Q] Z; // RE design matrix
  vector [K] Y[N]; // Response data
}

parameters {
  matrix [P, K-1] beta; 
  matrix [Q, K-1] eta; 
  vector<lower=0> [K-1] tau;
  real<lower=0> theta;
}

transformed parameters{
  matrix [N, K-1] lodds;
  matrix [N, K] odds;
  matrix [N, K] phi;
  matrix [Q, K-1] remat;
  
  for (k in 1:(K-1)){
    remat[, k] = eta[, k] * tau[k];
  }
  
  lodds = X*beta + Z*remat;
  for (n in 1:N){
    odds[n, 1:(K-1)] = exp(lodds[n, 1:(K-1)]);
    odds[n, K] = 1;
  }
  phi = odds/sum(odds);
}

model{
  // Priors
  for (p in 1:P){
    beta[p, ] ~ normal(0, 3);
  }
  
  tau ~ cauchy(0, 1);
  for (k in 1:(K-1)){
      eta[, k] ~ normal(0, 1);
  }
  
  theta ~ gamma(2, 0.1);
  
  // Likelihood
  for (n in 1:N){
    target += dirichlet_lpdf(Y[n] | (phi[n, ]*theta)');
  }
}

