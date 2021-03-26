data {
  int<lower=0> Np;
  int<lower=0> No;
  real<lower=0> lam1;
  real<lower=0> lam2;
  int<lower=0> pop;
  vector[Np] md;
  int hd[No];
}

transformed data{
  real md_o[No]; 
  for(i in 1:No){
    md_o[i] = sum(md[(((i-1)*52)+1):(i*52)]);
  }
  print(md_o);
}

parameters {
  vector<lower=0, upper=1>[No] zeta;
  real a;
  real b;
}

transformed parameters{
  vector<lower=0>[No] inf_m;
  for(i in 1:No){
    inf_m[i] = sum(md[(((i-1)*52)+1):(i*52)]) * zeta[i];
  }
  //print(inf_m);
  //print(zeta);
}

model {
  for(i in 1:No){
    hd[i] ~ poisson(pop*lam2*inf_m[i]/md_o[i]);
  }
  zeta ~ normal(a,b);
}

