derivative <- tsiR::derivative

## taken from the tsiR package
## need to be careful..
reconstruct <- function(cases,
						births,
						sigmamax=3){
	cumcases <- cumsum(cases)
	cumbirth <- cumsum(births)
	
	X <- cumcases
	Y <- cumbirth
	
	x <- seq(X[1], X[length(X)], length=length(X))
	y <- approxfun(X, Y)(x)
	y[1] <- y[2] - (y[3]-y[2])
	
	sigvec <- seq(sigmamax,0,-0.1)
	
	for(it in 1:length(sigvec)){
		Yhat <- predict(gausspr(x,y,variance.model=TRUE,fit=TRUE,tol=1e-7,
								var=9.999999999999999999e-3,
								kernel="rbfdot",
								kpar=list(sigma=sigvec[it])),X)
		
		
		if(sigvec[it] <= min(sigvec)){
			## use the loess then
			print('gaussian regressian failed -- switching to loess regression')
			Yhat <- predict(loess(y~x,se=T,family='gaussian',degree=1,model=T),X)
		}
		
		Z <- Y - Yhat
		rho <- derivative(X,Yhat)
		if(length(which(rho<=1))==0){
			break()
		}
	}
	
	adj.rho <- rho
	list(
		rho=adj.rho,
		Z=Z
	)
}

fitfun <- function(log.theta,
				   log.rho,
				   log.tau,
				   log.scale,
				   distmat,
				   Imat,
				   Smat,
				   popsize) {
	theta <- exp(log.theta)
	rho <- exp(log.rho)
	tau <- exp(log.tau)
	scale <- exp(log.scale)
	
	invdist <- 1/(distmat * scale)
	diag(invdist) <- 0
	
	## M[i, j] represents i to j movement
	M <- theta * t(popsize^tau * invdist^rho)
	
	diag(M) <- 1-rowSums(M)
	
	hatincmat <- Imat %*% M
	
	Inewmat <- tail(Imat, -1)
	
	Ipredmat <- head(hatincmat, -1)
	
	Spredmat <- head(Smat, -1)
	
	period <- rep(1:26, 100)[1:nrow(adjincmat)]
	
	fitdata <- data.frame(
		logInew=log(c(Inewmat)+1),
		logI=log(c(Ipredmat)+1),
		logS=log(c(Spredmat)),
		region=rep(names(popsize), each=nrow(Inewmat)),
		period=rep(head(period, -1), ncol(Inewmat)),
		logpop=rep(log(popsize), each=nrow(Inewmat))
	)
	
	fitdata$off <- fitdata$logS-fitdata$logpop
	
	gamm4fit <- try(gamm4(logInew ~ t2(period, region, k=26, bs=c("cc", "re"))  
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
				   popsize,
				   verbose=FALSE) {
	fit <- fitfun(log.theta, log.tau, log.rho, log.scale, distmat, Imat, Smat, popsize)
	
	gamm4fit <- fit$gamm4fit
	
	nll <- ifelse(inherits(gamm4fit, "try-error"), Inf, -logLik(gamm4fit$mer)[[1]])
	
	if (verbose) {
		print(paste("log.theta:", prettyNum(log.theta, digits=4), 
					"log.rho:", prettyNum(log.rho, digits=4), 
					"nll:", prettyNum(nll, digits=4)))
	}
	
	nll
}

