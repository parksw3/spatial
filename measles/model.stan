data {
	int ncity; // number of cities
	int nobs; // number of observations per city; this should be technically nobs - 1
	int nbasis;
	matrix[ncity, nobs] Inew; // tail(Imat, -1)
	matrix[ncity, nobs] Iprev; // head(Imat, -1)
	matrix[ncity, nobs] Zmat;
	vector[ncity] popsize; // population size
	matrix[nbasis, nobs] Bmat; // basis matrix
}

parameters {
	real logit_theta; // mixing proportion
	vector[ncity] logit_sprop; // proportion susceptible on a logit scale
	real mu_sprop;
	real<lower=0> sigma2_sprop;
	
	matrix[ncity, nbasis] amat;
	vector[nbasis] mu_a;
	real<lower=0> sigma2_a;
	
	vector<lower=0>[ncity] alpha;
	real<lower=0> mu_alpha;
	real<lower=0> sigma2_alpha;
	
	real phi;
}

transformed parameters {
	real<lower=0, upper=1> theta;
	vector<lower=0, upper=1>[ncity] sprop;
	vector[ncity] sbar;
	matrix[ncity, ncity] m;
	matrix[ncity, nobs] Ipred;
	matrix[ncity, nobs] logIhat;
	
	theta = inv_logit(logit_theta);
	sprop = inv_logit(logit_sprop);
	
	sbar = sprop * popsize;
	
	Smat = sbar + Zmat;
	
	Ipred = m * Iprev;
	
	for (i in 1:city) {
		logIhat[i,] = amat[i,] * Bmat + log(sbar[i] + Zmat[i,]) + alpha[i] * log(Ipred[i,]);
	}
}

model {
	logit_theta ~ normal(0, 20);
	
	logit_sprop ~ normal(mu_sprop, sqrt(sigma2_sprop));
	mu_sprop ~ normal(0, 20);
	sigma2_sprop ~ inv_gamma(1, 1);
	
	for (i in 1:nbasis) {
		amat[,i] ~ normal(mu_a[i], sqrt(sigma2_a));
	}
	
	mu_a ~ normal(0, 20);
	sigma2_a ~ inv_gamma(1, 1);
	
	alpha ~ normal(mu_alpha, sqrt(sigma2_alpha));
	mu_alpha ~ normal(1, 0.1);
	sigma2_alpha ~ inv_gamma(1,1);
	
	for (i in 1:city) {
		Inew[i,] ~ neg_binomial_2_log(logIhat[i,], phi);
	}
	phi ~ gamma(2, 1/30);
}
