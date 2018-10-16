data {
	int ncity; // number of cities
	int nobs; // number of observations per city; this should be technically nobs - 1
	int nbasis;
	int Inew[ncity, nobs]; // c(tail(Imat, -1))
	matrix[ncity, nobs] Iprev; // head(Imat, -1)
	matrix[ncity, nobs] Zmat;
	matrix[ncity, nobs] popmat; // population size head(popmat, -1)
	matrix[ncity, nobs] logpopmat; // population size head(popmat, -1)
	matrix[nbasis, nobs] Bmat; // basis matrix
	matrix[ncity, ncity] M;
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
	matrix[ncity, nobs] sbar;
	matrix[ncity, ncity] m; // m[i,j] models j->i
 	matrix[ncity, nobs] Ipred;
	matrix[ncity, nobs] logIhat;
	
	theta = inv_logit(logit_theta);
	sprop = inv_logit(logit_sprop);
	
	m = rep_matrix(theta, ncity, ncity) .* M;
	
	for (i in 1:ncity) {
		m[i,i] = 1 - sum(m[,i]);
	}
	
	Ipred = m * Iprev;
	
	for (i in 1:ncity) {
		sbar[i,] = rep_vector(sprop[i], nobs)' .* popmat[i,];
		
		logIhat[i,] = amat[i,] * Bmat + log(sbar[i,] + Zmat[i,]) + alpha[i] * log(Ipred[i,]) - logpopmat[i,];
	}
}

model {
	logit_theta ~ normal(-20, 5);
	
	logit_sprop ~ normal(mu_sprop, sqrt(sigma2_sprop));
	mu_sprop ~ normal(0, 10);
	sigma2_sprop ~ inv_gamma(1, 1);
	
	for (i in 1:nbasis) {
		amat[,i] ~ normal(mu_a[i], sqrt(sigma2_a));
	}
	
	mu_a ~ normal(log(34), 2);
	sigma2_a ~ inv_gamma(1, 1);
	
	alpha ~ normal(mu_alpha, sqrt(sigma2_alpha));
	mu_alpha ~ normal(1, 0.2);
	sigma2_alpha ~ inv_gamma(1,1);
	
	for (i in 1:ncity) {
		Inew[i,] ~ neg_binomial_2_log(logIhat[i,], phi);
	}
	
	phi ~ gamma(2, 0.025);
}
