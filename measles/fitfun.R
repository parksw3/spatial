fitfun <- function(log.theta,
				   log.rho,
				   log.tau,
				   log.scale,
				   distmat,
				   Imat,
				   Smat,
				   popmat) {
	theta <- exp(log.theta)
	rho <- exp(log.rho)
	tau <- exp(log.tau)
	scale <- exp(log.scale)
	
	invdist <- 1/(distmat * scale)
	diag(invdist) <- 0
	
	## M[i, j] represents i to j movement
	M <- theta * t(popmat[1,]^tau * invdist^rho)
	
	diag(M) <- 1-rowSums(M)
	
	hatincmat <- Imat %*% M
	
	Inewmat <- tail(Imat, -1)
	
	Ipredmat <- head(hatincmat, -1)
	
	Spredmat <- head(Smat, -1)
	
	period <- rep(1:52, 100)[1:nrow(Inewmat)]
	
	fitdata <- data.frame(
		logInew=log(c(Inewmat)+1),
		logI=log(c(Ipredmat)),
		logS=log(c(Spredmat)),
		region=rep(colnames(Inewmat), each=nrow(Inewmat)),
		period=rep(period, ncol(Inewmat)),
		logpop=c(head(log(popmat), -1))
	)
	
	fitdata$off <- fitdata$logS-fitdata$logpop
	
	gamm4fit <- try(gamm4(logInew ~ t2(period, region, k=52, bs=c("cc", "re"))  
						 + logI + offset(off),
						 random=~(logI||region),
						 data=fitdata))
	
	list(
		gamm4fit=gamm4fit,
		theta=theta,
		rho=rho,
		scale=scale,
		M=M,
		distmat=distmat,
		Imat=Imat,
		Smat=Smat
	)
}

nllfun <- function(log.theta,
				   log.rho,
				   log.tau,
				   log.scale,
				   distmat,
				   Imat,
				   Smat,
				   popmat,
				   verbose=FALSE) {
	fit <- fitfun(log.theta, log.tau, log.rho, log.scale, distmat, Imat, Smat, popmat)
	
	gamm4fit <- fit$gamm4fit
	
	nll <- ifelse(inherits(gamm4fit, "try-error"), Inf, -logLik(gamm4fit$mer)[[1]])
	
	if (verbose) {
		print(paste("log.theta:", prettyNum(log.theta, digits=4), 
					"log.rho:", prettyNum(log.rho, digits=4), 
					"nll:", prettyNum(nll, digits=4)))
	}
	
	nll
}
