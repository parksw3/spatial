simulate.sir <- function(betamat, ## transmission matrix ncity by nperiod
					 alpha,
					 mixmat, ## m[i, j] represents j -> i
					 I0, S0,
					 nsim=10,
					 method=c("nbinom", "deterministic", "poisson"),
					 popmat,
					 birthmat, 
					 rhomat, ## rhomat needs to be a column longer
					 phi) { 
	method <- match.arg(method)
	
	pp <- ncol(betamat)
	
	ncity <- nrow(popmat)
	
	mfun <- switch(method,
				   deterministic=function(x, y) x,
				   poisson=function(x, y) rpois(length(x), x),
				   nbinom=function(x, y) rnbinom(length(x), mu=x, size=y))
	
	rfun <- switch(method,
				   function(x, p) rbinom(length(x), round(x), p),
				   deterministic=function(x,p) x*p) 
	
	pfun <- switch(method,
				   function(m,I) {
				   	if (ncity > 1) {
				   		rowSums(sapply(1:ncity, function(x){rmultinom(1, I[x], m[,x])}))
				   	} else {
				   		I
				   	}
				   },
				   deterministic=function(m,I) m %*% I)
	
	bfun <- switch(method,
				   function(x) rpois(length(x), x),
				   deterministic=function(x) x
	)
	
	ncity <- nrow(betamat)
	nperiod <- ncol(betamat)
	
	reslist <- vector('list', nsim)
	
	for (i in 1:nsim) {
		simmatI <- matrix(0, nrow=nrow(popmat), ncol=ncol(popmat)+1)
		simmatC <- matrix(0, nrow=nrow(popmat), ncol=ncol(popmat)+1)
		simmatS <- matrix(0, nrow=nrow(popmat), ncol=ncol(popmat)+1)
		
		simmatI[,1] <- I0
		simmatS[,1] <- S0
		simmatC[,1] <- rfun(I0, 1/rhomat[,1])
		
		for (t in 2:ncol(simmatI)) {
			period <- t %% pp
			
			if (period==0) period <- pp
			
			Iprev <- Iprev2 <- pfun(mixmat, simmatI[,t-1])
			
			Iprev2[Iprev2==0] <- 1
			
			if (missing(phi)) {
				size <- Iprev2
			} else {
				size <- phi
			}
			
			if (all(Iprev==0)) {
				for (j in t:ncol(simmatI)) {
					simmatS[,j] <- simmatS[,t-1]
				}
				
				break
			}
			
			incidence <- mfun(betamat[,period]*simmatS[,t-1]*(Iprev)^alpha/popmat[,t-1], size)
			
			ii <- pmin(simmatS[,t-1], incidence)
			
			simmatI[,t] <- ii
			simmatS[,t] <- simmatS[,t-1] - ii + bfun(birthmat[,t-1])
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
