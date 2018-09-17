runltsir <- function(data,
					 g_t,
					 IP=2,
					 sigmamax=3) {
	cumcases <- cumsum(data$cases)
	cumbirth <- cumsum(data$births)
	
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
		
		Z <- residual.cases(Yhat,Y)
		rho <- derivative(X,Yhat)
		if(length(which(rho<=1))==0){
			break()
		}
	}
	
	adj.rho <- rho
	
	period <- rep(1:(52/IP), round(nrow(data)+1))[1:(nrow(data)-2)]
	
	Iadjusted <- round(data$cases * adj.rho)
	
	lag_length <- length(g_t)
	
	Inew <- tail(Iadjusted,-lag_length) ## t+1
	
	lagImat <- matrix(0, nrow=length(Inew), ncol=lag_length)
	
	for (i in 1:lag_length) {
		lagImat[,i] <- 
			Iadjusted[1:length(Inew) + (i-1)]
	}
	
	Zminus <- Z[1:length(Inew) + (i-1)]
	
	g_mat <- matrix(0, lag_length, lag_length)
	
	diag(g_mat) <- rev(g_t_adj)
	
	lagImat_adj <-lagImat %*% g_mat
	
	pop <- data$pop
	minSmean <- max(0.01*pop,-(min(Z)+1))
	Smean <- exp(seq(log(minSmean), log(mean(pop)), length=500))
	nll <- rep(NA, 500)
	
	nrow <- length(Inew)
	ncol <- 26
	
	for (i in 1:500) {
		Smean0 <- Smean[i]
		
		mm <- matrix(0, nrow=ncol, ncol=nrow)
		
		for (j in 1:lag_length) {
			diag <- (1:nrow+(j-1))%%ncol 
			diag[diag==0] <- 26
			
			mm[0:(nrow-1)*26+diag] <- lagImat_adj[,j] * (Zminus + Smean0)
		}
		
		mm_adj <- t(mm)
		
		fit0 <- try(glm.fit(mm_adj, Inew, family=Gamma("identity")))
		
		nll[i] <- -try(stats:::logLik.lm(fit0))
	}
	
	k <- which.min(nll)
	
	Smean0 <- Smean[k]
	
	mm <- matrix(0, nrow=ncol, ncol=nrow)
	
	for (j in 1:lag_length) {
		diag <- (1:nrow+(j-1))%%ncol 
		diag[diag==0] <- 26
		
		mm[0:(nrow-1)*26+diag] <- (Zminus + Smean0) * lagImat_adj[,j] 
	}
	
	mm_adj <- t(mm)
	
	fit <- glm.fit(mm_adj, Inew, family=Gamma("identity"))
}

runltsir_NB <- function(data,
						IP=2,
						sigmamax=3,
						spline=FALSE,
						predict=FALSE,
						nsim=10,
						method="pois",
						add.noise.sd = 0, mul.noise.sd = 0) {
	cumcases <- cumsum(data$cases)
	cumbirth <- cumsum(data$births)
	
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
		
		Z <- residual.cases(Yhat,Y)
		rho <- derivative(X,Yhat)
		if(length(which(rho<=1))==0){
			break()
		}
	}
	
	adj.rho <- rho
	
	period <- rep(1:(52/IP), round(nrow(data)+1))[2:(nrow(data)-1)]
	
	Iadjusted <- round(data$cases * adj.rho)
	
	lag_length <- 2
	
	Inew <- tail(Iadjusted,-lag_length) ## t+1
	
	lagImat <- matrix(0, nrow=length(Inew), ncol=lag_length)
	
	for (i in 1:lag_length) {
		lagImat[,i] <- 
			Iadjusted[1:length(Inew) + (i-1)]
	}
	
	Zminus <- Z[1:length(Inew) + (i-1)]
	
	oo <- optim(
		c(sbar=log(2*abs(min(Zminus))), epsilon=qlogis(0.1)),
		runltsir_NB_obj, 
		Inew=Inew,
		lagImat=lagImat,
		Zminus=Zminus,
		period=period,
		spline=spline
	)
	
	sbar <- exp(oo$par[1])
	epsilon <- plogis(oo$par[2])
	
	g_mat <- matrix(0, 2, 2)
	
	diag(g_mat) <- c(epsilon, 1-epsilon)
	
	lagImat_adj <- lagImat %*% g_mat
	
	lSpred <- log(sbar + Zminus)
	
	lIpred <- log(rowSums(lagImat_adj))
	
	fitdata <- data.frame(
		period=period,
		lSpred=lSpred,
		lIpred=lIpred
	)
	
	if (spline) {
		fit <- gam(Inew~ s(period, bs="cc") + offset(lSpred) + lIpred,
				   family=nb,
				   data=fitdata)
		
		beta <- exp(predict(
			fit, 
			newdata=data.frame(period=1:26, lIpred=0, lSpred=0)
		))
		alpha <- coef(fit)[2]
	} else {
		fit <- MASS::glm.nb(Inew~ -1 + as.factor(period) + offset(lSpred) + lIpred,
							data=fitdata)
		beta <- exp(coef(fit)[1:26])
		alpha <- coef(fit)[27]
	}
	
	if (predict) {
		S <- rep(0,length(data$cases))
		I <- rep(0,length(data$cases))
		
		S_start <- sbar + Z[1:2]
		
		head(lagImat_adj)
		
		res <- matrix(0,length(data$cases),nsim)
		Sres <- matrix(0,length(data$cases),nsim)
		
		for(ct in 1:nsim){
			
			S <- rep(0,length(data$cases))
			I <- rep(0,length(data$cases))
			
			S[1:2] <- S_start[1:2]
			I[1:2] <- adj.rho[1:2]*data$cases[1:2]
			
			for (t in 3:(nrow(data))){
				
				lambda <- min(S[t-1], 
							  unname(beta[period[t-2]] * S[t-1] * (I[t-1] * (1-epsilon) + I[t-2] * epsilon)^alpha))
				
				#if(lambda < 1 || is.nan(lambda) == T){lambda <- 0}
				if(is.nan(lambda) == T){lambda <- 0}
				
				if(method == 'deterministic'){
					I[t] <- lambda * rnorm( n = 1, mean = 1, sd=mul.noise.sd)
					if(I[t] < 0 && lambda >= 0 ){
						warning('infected overflow  -- reduce multiplicative noise sd')
					}
				}
				if(method == 'negbin'){
					I[t] <- rnbinom(n=1,mu=lambda,size=I[t-1]+1e-10)
				}
				if(method == 'pois'){
					I[t] <- rpois(n=1,lambda=lambda)
				}
				
				S[t] <- max(S[t-1] + data$births[t-1] - I[t] + rnorm(n=1,mean=0,sd=add.noise.sd),0)
				
				if(S[t] < 0 && (S[t-1] + data$births[t-1] - I[t]) >0 ){
					warning('susceptible overflow  -- reduce additive noise sd')
				}
			}
			res[,ct] <- I / adj.rho
			Sres[,ct] <- S
		}
		
	} else {
		res <- NULL
	}
	
	list(
		fit=fit,
		epsilon=epsilon,
		sbar=sbar,
		res=res,
		Sres=Sres
	)
}

runltsir_NB_obj <- function(par,
							Inew,
							lagImat,
							Zminus,
							period,
							spline) {
	sbar <- exp(par[1]); epsilon <- plogis(par[2])
	
	g_mat <- matrix(0, 2, 2)
	
	diag(g_mat) <- c(epsilon, 1-epsilon)
	
	lagImat_adj <- lagImat %*% g_mat
	
	Ipred <- rowSums(lagImat_adj)
	
	if (spline) {
		fit <- gam(Inew~ s(period, bs="cc") + offset(log(sbar + Zminus)) + log(Ipred),
				   family=nb)
	} else {
		fit <- MASS::glm.nb(Inew~ -1 + as.factor(period) + offset(log(sbar + Zminus)) + log(Ipred))
	}

	-logLik(fit)[[1]]
}

