functions{
  real logLikeBDcoalTimes_lpdf(real[] t, real lambda, real mu, real lgRho){
    int numCoal;
    real ll;
    numCoal=size(t);

    // Define the log likelihood
    ll=0.0;
    ll = ll + log(numCoal + 1); // First term in place of doing a factorial (I dont know how)
    ll = ll + log(lambda - mu) + (numCoal)*(log(lambda) + lgRho) - (lambda - mu)*max(t);
    ll = ll - log(exp(lgRho)*lambda + (lambda*(1-exp(lgRho))-mu)*exp(-(lambda-mu)*max(t)));
    for(k in 1:numCoal){
      ll = ll + log(k); // Subsequent terms of factorial
      ll = ll + 2*log(lambda-mu) - (lambda-mu)*t[k];
      ll = ll - 2*log(exp(lgRho)*lambda + (lambda*(1-exp(lgRho)) - mu)*exp(-(lambda-mu)*t[k]));
    }

    return(ll);
  }
}

data{
  int<lower=1> nCoal; // Number of coalescence times
  real t[nCoal];  // Coalescence times
  real upperLambda; // Max growth rate allowed
}

parameters {
  real<lower=0, upper=upperLambda> lambda;
  real<lower=0, upper=lambda> mu;
  real<lower=-1000, upper=0> lgRho;
}

model { // Lack of specified priors indicates a uniform prior on the parameter bounds
  t ~ logLikeBDcoalTimes(lambda, mu, lgRho);
}
