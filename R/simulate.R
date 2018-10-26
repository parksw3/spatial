simulate_sir <- function(alpha,
						 beta,
						 tmax,
						 birth,
						 I0, S0, N,
						 sim.method=c("deterministic", "poisson", "nbinom")){
	period <- length(beta)
	
	sim.method <- match.arg(sim.method)
	
	mfun <- switch(sim.method,
				   deterministic=function(x, y) x,
				   poisson=function(x, y) rpois(length(x), x),
				   nbinom=function(x, y) rnbinom(length(x), mu=x, size=y))
	
	Ivec <- Svec <- rep(0, tmax)
	
	Ivec[1] <- I0
	Svec[1] <- S0
	
	for (i in 2:tmax) {
		j <- i %% period
		
		if (j==0) j <- period
		
		bb <- beta[j]
		
		inc <- min(Svec[i-1], bb * Ivec[i-1]^alpha * Svec[i-1]/N)
		
		Ivec[i] <- inc
		Svec[i] <- Svec[i-1] - inc + birth[i-1]
		
	}
	
	data.frame(
		I=Ivec,
		S=Svec
	)
	
}
