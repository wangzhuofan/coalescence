data {
  int<lower = 0> n;
  int<lower = 0> px;
  int<lower = 0> py;
  
  simplex[px] x[n];
  simplex[py] y[n];
}

parameters {
  real<lower=0,upper=1> c;
  vector<lower=0>[px] a;
  vector<lower=0>[py] b;
  
  simplex[py] m[px];
}

model {
  matrix[n,px] xm;
  matrix[n,py] ym;
  matrix[n,py] bound;
  matrix[n,py] rr;
  matrix[px,py] mm;
  for(i in 1:px) mm[px] = m[i]';
  for(i in 1:n) {xm[i] = x[i]';ym[i] = y[i]';}
  bound = xm*mm;
  // for(i in 1:n) {
  //   for(j in 1:py) {
  //     bound[i,j]=0;
  //     for(k in 1:px) {
  //       bound[i,j]=bound[i,j]+x[i,k]*m[k,j];
  //     }
  //     bound[i,j]=y[i,j]/bound[i,j];
  //   }
  // }
  bound = ym ./bound;
  c~beta(1,1);
  a~gamma(1,1);
  b~gamma(1,1);
  for(j in 1:px) m[j]~dirichlet(rep_vector(1,py));
  rr = ym-c*min(bound)*xm*mm;
  for(i in 1 : n) {
    x[i] ~ dirichlet(a);
    // rr[i] = y[i]';
    // for(j in 1:py) {
    //   for(k in 1:px) {
    //     rr[i,j]=rr[i,j]-c*min(bound)*x[i,k]*m[k,j];
    //   }
    // }
    rr[i]=rr[i]/(1-c*min(bound));
    
    rr[i]'~dirichlet(b);
    target+= 0-py*log(1-c*min(bound));
  }
}
