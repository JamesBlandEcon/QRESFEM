functions {
  
  
}



data {
  int<lower=0> N;
  int<lower=0> nstrg;
  vector[nstrg] follow[N];
  int treatmentID[N];
  int<lower=0> nTreatment;
  vector[N] count;
  
  vector[nstrg] prior_MIX;
  
  vector[N] UseData;
  
  
}
transformed data {
  
  matrix[N,nstrg] COUNT = rep_matrix(count,nstrg);
  
}


parameters {
  
  array[nTreatment] simplex[nstrg] MIX;
  vector<lower=0,upper=0.5>[nTreatment] epsilon;
  
  
  
}

transformed parameters {
    
    vector[N] log_lik;
  
   for (ii in 1:N) {
     
     int t = treatmentID[ii];
     
     log_lik[ii] = log_sum_exp(
       log(MIX[t])
                +follow[ii]*log(1.0-epsilon[t])+(count[ii]-follow[ii])*log(epsilon[t])
                );
      
    }
   
}


model {
    
    target += dirichlet_lpdf(MIX | prior_MIX);
    
    target += UseData.*log_lik;
    
   


  
}

generated quantities {

}
