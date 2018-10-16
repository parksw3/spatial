simfun <- function(tmat, ## transmission matrix ncity by nperiod
				   alpha,
				   mixmat, ## m[i, j] represents j -> i
				   I0, S0,
				   nsim=10,
				   method=c("nbinom", "deterministic", "poisson"),
				   popmat,
				   birthmat, ## rhomat needs to be a column longer
				   rhomat) { 
	method <- match.arg(method)
	
	mfun <- switch(method,
				   deterministic=function(x, y) x,
				   poisson=function(x, y) rpois(length(x), x),
				   nbinom=function(x, y) rnbinom(length(x), mu=x, size=y))
	
	rfun <- switch(method,
				   function(x, p) rbinom(length(x), round(x), p),
				   deterministic=function(x,p) x*p) 
	
	ncity <- nrow(tmat)
	nperiod <- ncol(tmat)
	
	reslist <- vector('list', nsim)
	
	for (i in 1:nsim) {
		simmatI <- matrix(0, nrow=nrow(popmat), ncol=ncol(popmat)+1)
		simmatC <- matrix(0, nrow=nrow(popmat), ncol=ncol(popmat)+1)
		simmatS <- matrix(0, nrow=nrow(popmat), ncol=ncol(popmat)+1)
			
		simmatI[,1] <- I0
		simmatS[,1] <- S0
		simmatC[,1] <- rfun(I0, 1/rhomat[,1])
		
		for (t in 2:ncol(simmatI)) {
			period <- t %% 52
			
			if (period==0) period <- 52
			
			Iprev <- mixmat %*% simmatI[,t-1]
			
			if (all(Iprev==0)) break
			
			incidence <- mfun(tmat[,period]*simmatS[,t-1]*(Iprev)^alpha/popmat[,t-1], Iprev)
			
			ii <- pmin(simmatS[,t-1], incidence)
			
			simmatI[,t] <- ii
			simmatS[,t] <- simmatS[,t-1] - ii + birthmat[,t-1]
			simmatC[,t] <- rfun(ii, 1/rhomat[,t])
		}
		
		reslist[[i]] <- list(
			I=simmatI,
			S=simmatS,
			C=simmatC
		)
	}
	reslist
}
