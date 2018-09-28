library(kernlab)
library(dplyr)
library(gamm4)

measles_data <- read.csv("measlesUKUS.csv")

derivative <- tsiR::derivative

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

measles_UK <- measles_data %>% 
	filter(country=="UK") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases))

measles_loc <- unique(measles_UK$loc)

measles_list <- measles_UK %>%
	split(as.character(.$loc))

sublist <- measles_list[c("LONDON", "BIRMINGHAM", "LIVERPOOL", "MANCHESTER", "BRISTOL")]

reconstruct_list <- sublist %>%
	lapply(function(data) reconstruct(data$cases, data$rec))

Zmat <- reconstruct_list %>%
	lapply("[[", "Z") %>%
	do.call(what="cbind")

rhomat <- reconstruct_list %>%
	lapply("[[", "rho") %>%
	do.call(what="cbind")

popsize <- sublist %>%
	lapply("[[", "pop") %>%
	do.call(what="cbind") %>%
	"["(1,)

## assumption
sbar <- 0.035 * popsize

Smat <- Zmat + matrix(rep(sbar, nrow(Zmat)), ncol=ncol(Zmat), byrow=TRUE)

any(Smat < 0)

incmat <- sublist %>%
	lapply("[[", "cases") %>%
	do.call(what="cbind")

## true (?) incidence
adjincmat <- round(incmat * rhomat)

period <- rep(1:26, 100)[1:nrow(adjincmat)]

distmat <- sublist %>%
	lapply(filter, decimalYear==min(decimalYear)) %>%
	bind_rows(.id="city") %>%
	select(lon, lat) %>%
	dist(diag=TRUE, upper=TRUE) %>%
	as.matrix



fitfun <- function(theta,
				   tau) {
	## M[i, j] represents i to j movement
	
	M <- theta * t(popsize^tau * distmat)
	
	diag(M) <- 1-rowSums(M)
	
	hatincmat <- adjincmat %*% M
	
	Inewmat <- tail(adjincmat, -1)
	
	Ipredmat <- head(hatincmat, -1)
	
	Spredmat <- head(Smat, -1)
	
	fitdata <- data.frame(
		logInew=log(c(Inewmat)+1),
		logI=log(c(Ipredmat)+1),
		logS=log(c(Spredmat)),
		region=rep(names(popsize), each=nrow(Inewmat)),
		period=rep(head(period, -1), ncol(Inewmat)),
		logpop=rep(log(popsize), each=nrow(Inewmat))
	)
	
	fitdata$off <- fitdata$logS-fitdata$logpop
	
	lmerfit <- gamm4(logInew ~ t2(period, region, k=26, bs=c("cc", "re")) 
					 + logI + offset(off),
					 data=fitdata)
	
}
