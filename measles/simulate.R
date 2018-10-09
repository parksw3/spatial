simfun <- function(fitted,
				   nsim=10,
				   method=c("deterministic", "poisson", "nbinom"),
				   popmat,
				   birthmat) {
	method <- match.arg(method)
	
	mfun <- switch(method,
				   deterministic=function(x) x,
				   poisson=function(x) rpois(length(x), x),
				   nbinom=function(x) rnbinom(length(x), mu=x, size=x))
	
	pp <- exp(predict(fitted$gamm4fit$gam,
				  newdata=data.frame(
				  	period=rep(1:52, ncol(fitted$distmat)),
				  	logI=0,
				  	off=0,
				  	region=rep(colnames(fitted$Imat), each=52) 
				  )))
	
	alpha <- fixef(fitted$gamm4fit$mer)[2] + ranef(fitted$gamm4fit$mer)$region$logI
	
	trans <- matrix(pp, nrow=52)
	
	reslist <- vector('list', nsim)
	
	mixmat <- fitted$M
	
	for (i in 1:nsim) {
		simmatI <- matrix(0, nrow=nrow(fitted$Imat), ncol=ncol(fitted$Imat))
		simmatS <- matrix(0, nrow=nrow(fitted$Imat), ncol=ncol(fitted$Imat))
			
		simmatI[1,] <- fitted$Imat[1,]
		simmatS[1,] <- fitted$Smat[1,]
		
		for (t in 2:nrow(simmatI)) {
			period <- t %% 52
			
			if (period==0) period <- 52
			
			incidence <- mfun(trans[period,]*simmatS[t-1,]*(simmatI[t-1,] %*% mixmat)^alpha/popmat[t-1,])
			
			ii <- pmin(simmatS[t-1,], incidence)
			
			simmatI[t,] <- ii
			simmatS[t,] <- simmatS[t-1,] - ii + birthmat[t-1,]
		}
		
		reslist[[i]] <- list(
			I=simmatI,
			S=simmatS
		)
	}
	
}
