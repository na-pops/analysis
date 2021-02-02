/* Where I left off:
  -------------------
  Trying to figure out how to index by species, then by sampling event (or sample
  event X time band). Basically, I want to grab all counts belonging to species
  x, which is easy enough because the original data are ordered by species, so
  it is just a matter of using the segment() function for each species. But where
  I am having a brain fart is that there are j_i time bands per sample event i,
  where length of j_i are not all equal. Anyway, this is tomorrow, or Monday's 
  issue.

*/


data {
  int<lower = 1> n_samples;           // total number of sampling events i
  int<lower = 1> n_covariates;        // total number of covariates
  int<lower = 1> n_species;           // total number of species being modelled
  int<lower = 2> max_intervals;       // maximum number of intervals being considered
  int samples_per_species[n_species]; // number of samples per species
  int abund_per_band[n_samples, max_intervals];// abundance in time band j for sample i
  vector[n_samples] abund_per_sample; // total abundnace for sample i
  int bands_per_sample[n_samples]; // number of time bands for sample i
  matrix[n_samples, max_intervals] max_time; // max time duration for time band j
  matrix[n_samples, n_covariates] X;  // matrix of covariates
}

parameters {
  row_vector[n_covariates] mu;           // mean vector
  vector<lower = 0>[n_covariates] tau;   // scale vector for Sigma
  corr_matrix[n_covariates] Omega;       // correlation matrix for Sigma
  matrix[n_species, n_covariates] gamma; // coefficients of interest
}

model {
  vector[n_samples] log_phi;             // singing rate
  matrix[n_samples, max_intervals] Pi;   // probabilities
  
  Omega ~ lkj_corr(2);
  tau ~ cauchy(0, 2.5);
  mu ~ normal(0, 0.5);
  
  Pi = rep_matrix(0, n_samples, max_intervals);
  
  for (s in 1:n_species)
  {
    gamma[,s] ~ multi_normal(mu, quad_form_diag(Omega, tau));
  }
  log_phi = rows_dot_product(X, to_matrix(gamma));
  
  for (i in 1:n_samples)
  {
    for (j in 2:bands_per_sample[i])
    {
      Pi[i,j] = (exp(-max_time[i,j-1] * exp(log_phi[i])) - 
                 exp(-max_time[i,j] * exp(log_phi[i]))) / 
                (1 - exp(-max_time[i,bands_per_sample[i]] * exp(log_phi[i])));
    }
    Pi[i,1] = 1 - sum(Pi[i,]);
    
    abund_per_band[i,] ~ multinomial(to_vector(Pi[i,]));
  }

}

