derivative <- tsiR::derivative

reconstruct_scam <- function(cases,
						births,
						k=10,
						m=6){
	cumcases <- cumsum(cases)
	cumbirth <- cumsum(births)
	
	dd <- data.frame(
		x=cumcases,
		y=cumbirth
	)
	
	fit <- scam(y~s(x, k=k, bs="mpi", m=m), data=dd)
	
	adj.rho <- derivative(dd$x, predict(fit))
	
	list(
		rho=adj.rho,
		Z=fit$residuals,
		Yhat=predict(fit)
	)
}

reconstruct_gauss <- function(cases,
							  births,
							  sigmamax=3) {
	cumcases <- cumsum(cases)
	cumbirth <- cumsum(births)
	
	x <- cumcases
	y <- cumbirth
	
	
	sigvec <- seq(sigmamax,0,-0.1)
	for(it in 1:length(sigvec)){
		
		Yhat <- predict(gausspr(x,y,variance.model=T,fit=T,tol=1e-7,
								var=9.999999999999999999e-3,
								kernel="rbfdot",
								kpar=list(sigma=sigvec[it])),x)
		
		
		if(sigvec[it] <= min(sigvec)){
			## use the loess then
			print('gaussian regressian failed -- switching to loess regression')
			Yhat <- predict(loess(y~x,se=T,family='gaussian',degree=1,model=T),x)
		}
		
		Z <- y-Yhat
		rho <- derivative(x,Yhat)
		if(length(which(rho<=1))==0){
			break()
		}
	}
	
	list(
		rho=rho,
		Z=Z,
		Yhat=Yhat
	)
}
