functions {
  real VST(real g,real l,real delta,real epsilon,matrix numFollow, int[,] nextState,int firstState) {
    
    vector[4] V;
    matrix[4,4] h = exp(numFollow*log(1.0-epsilon)+(2.0-numFollow)*log(epsilon));
    matrix[4,4] H = rep_matrix(0,4,4);
    
    vector[4] u = to_vector({0,g+1.0,-l,1.0});
    for (ss in 1:4) {
      for (oo in 1:4) {
        H[ss,nextState[ss,oo]] = H[ss,nextState[ss,oo]]+h[ss,oo];
      }
    }
    
    V = (1-delta)*(diag_matrix(rep_vector(1,4))-delta*H)\u;
    
    return V[firstState];
    
    
  }
  
  vector solveQRE(matrix U,matrix A,vector lgrid, int nc,data real ftol) {
    int sizeLgrid = num_elements(lgrid);
    int nStrg = dims(U)[1];
    
    vector[nStrg] q0 = rep_vector(log(1.0/nStrg),nStrg);
    
    for (ll in 2:sizeLgrid) {
      
      real lm = 0.5*(lgrid[ll]-lgrid[ll-1]);
      vector[nStrg] Hlambda = append_row(-A*U*exp(q0),0);
      matrix[nStrg,nStrg] Hq = append_row(A-lm*A*U*diag_matrix(exp(q0)),
                                          -rep_row_vector(1,nStrg)*diag_matrix(exp(q0)));
      
      q0 = q0-Hq\Hlambda * (lgrid[ll]-lgrid[ll-1]);
      
      
      for (cc in 1:nc) {
        vector[nStrg] H;
        
        H  = append_row(A*q0-lgrid[ll]*A*U*exp(q0),1.0-sum(exp(q0)));
        if (sqrt(H'*H)>ftol) {
        Hq = append_row(A-lgrid[ll]*A*U*diag_matrix(exp(q0)),
                        -rep_row_vector(1,nStrg)*diag_matrix(exp(q0)));
        
        q0 = q0-Hq\H;
        } else {
          break;
        }
        
        if (cc==nc) {
          print("maxiter reached!");
        }
      }
      
      
      
    }
    
    return q0;
  }
  
}



data {
  int<lower=0> N;
  int<lower=0> nstrg;
  vector[nstrg] follow[N];
  int treatmentID[N];
  int<lower=0> nTreatment;
  
  real g[nTreatment]; 
  real l[nTreatment];
  real delta[nTreatment];
  
  vector[N] count;
  matrix[4,4] numFollow[nstrg,nstrg];
  int nextState[nstrg,nstrg,4,4];
  int firstState[nstrg,nstrg];
  
  real prior_lambda[2];
  
  int lGridSize; 
  int nc;
  
  real ftol;
  vector[N] UseData;
  
  
}
transformed data {
  matrix[nstrg-1,nstrg] A;
  for (rr in 1:(nstrg-1)) {
    for (cc in 1:nstrg) {
      A[rr,cc] = 0;
      if (cc==(rr+1)) {
        A[rr,cc]=-1;
      }
      if (cc==1) {
        A[rr,cc]=1;
      }
    }
  }
  
  print("This is QRE-SFEM!");
  print("James R. Bland tried this in 2025 and it didin't work particularly well when applied to Dal Bo and Frechette (2011). As such, you may want to think twice about implementing it for yourself. Go on, put it back in the box and just estimate the classic SFEM!");
  
}


parameters {
  
  real<lower=0> lambda;
  real<lower=0,upper=0.5> epsilon;
  
}

transformed parameters {
   vector[nstrg] logmix[nTreatment];
  
  for (tt in 1:nTreatment) {
    matrix[nstrg,nstrg] U;
    vector[lGridSize] lgrid;
    for (rr in 1:nstrg) {
      for (cc in 1:nstrg) {
        U[rr,cc] = VST(g[tt],l[tt],delta[tt],epsilon,numFollow[rr,cc], nextState[rr,cc,,],firstState[rr,cc]);
      }
    }
    for (ll in 1:lGridSize) {
      lgrid[ll] =lambda*(1.0/(1.0-0.5*(ll-1.0)/(lGridSize-1.0))-1);
    }
    logmix[tt] = solveQRE(U,A,lgrid,nc,ftol);
    
   
    
  }
}


model {
  

  target += lognormal_lpdf(lambda | prior_lambda[1],prior_lambda[2]);
  
  
  for (ii in 1:N) {
    target += UseData[ii]*log_sum_exp((follow[ii]*log(1.0-epsilon)+(count[ii]-follow[ii])*log(epsilon))+logmix[treatmentID[ii]]);
  }
}

generated quantities {
  vector[N] log_lik;
  for (ii in 1:N) {
    log_lik[ii]=log_sum_exp((follow[ii]*log(1.0-epsilon)+(count[ii]-follow[ii])*log(epsilon))+logmix[treatmentID[ii]]);
  }
  
  vector[nstrg] MIX[nTreatment];
  
  for (tt in 1:nTreatment) {
    MIX[tt] = exp(logmix[tt]);
  }

}
