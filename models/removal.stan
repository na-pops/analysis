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
  int<lower = 1> n_counts;            // total number of counts (i x j ∀i)
  int<lower = 1> n_samples;           // total number of sampling events i
  int<lower = 1> n_covariates;        // total number of covariates
  int<lower = 1> n_species;           // total number of species being modelled
  
  vector[n_counts] abund_per_band;    // abundance in time band j for sample i
  vector[n_samples] abund_per_sample; // abundance for sample i, species s
  vector[n_samples] bands_per_sample; // number of time bands for sample i
  
  vector[n_species] samples_per_species; // number of samples per species
  vector[n_species] count_per_species;   // number of counts per species
  vector[n_samples] species;             // species for sample i
  
  vector[n_counts] max_time;          // max time duration for time band j
  
  matrix[n_samples, n_covariates] X;  // matrix of covariates
}

parameters {
  vector[n_covariates] mu;               // mean vector
  vector<lower = 0>[n_covariates] tau;   // scale vector for Sigma
  corr_matrix[n_covariates] Omega;       // correlation matrix for Sigma
  matrix[n_species, n_covariates] gamma; // coefficients of interest
  vector[n_samples] log_phi;             // singing rate
  vector[n_counts] Pi;
}

model {
  Omega ~ lkj_corr(2);
  tau ~ cauchy(0, 2.5);
  mu ~ normal(0, 0.5);
  
  // Adapted from Stan User's Guide Ch 8.2 Ragged data structures
  int min_i = 1;
  int max_i = 0;
  int pos = 1;
  
  for (s in 1:n_species)
  {
    gamma[s,] ~ multi_normal(mu, quad_form_diag(Omega, tau));
    
    max_i = max_i + samples_per_species[s];
    for (i in min_i:max_i)
    {
      log_phi[i,s] = X[i,] * gamma[s,]';
      
      
      sample_max_time = segment(max_time, pos, bands_per_sample[i]);
      for (j in 2:bands_per_sample[i])
      {
        
      }
      
      pos = pos + bands_per_sample[i];
    }
    
    min_i = min_i + samples_per_species[s];
  }

}

