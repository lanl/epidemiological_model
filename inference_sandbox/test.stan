data {
  int<lower=0> Np;
  int<lower=0> No;
  real<lower=0> lam1;
  real<lower=0> lam2;
  int<lower=0> pop;
  vector[Np] md;
  vector[Np] bd;
  int hd[No];
}

transformed data{
  real md_o[No]; 
  real bd_o[No];
  for(i in 1:No){
    md_o[i] = sum(md[(((i-1)*52)+1):(i*52)]);
    bd_o[i] = sum(bd[(((i-1)*52)+1):(i*52)]);
  }
  //print(md_o);
}

parameters {
  vector<lower=0, upper=1>[No] zeta;
  real a;
  real b;
}

transformed parameters{
  vector<lower=0>[No] inf_mo;
  for(i in 1:No){
    inf_mo[i] = md_o[i] * bd_o[i] * lam1 * zeta[i];
    if(inf_mo[i]>md_o[i])inf_mo[i] = md_o[i];
  }
}

model {
  for(i in 1:No){
    hd[i] ~ poisson(pop*lam2*inf_mo[i]);
  }
  zeta ~ normal(a,b);
}

