library(mgcv)
library(dplyr)
library(scam)
library(rstan)
source("../R/reconstruct.R")
source("fitfun.R")

measles_data <- read.csv("measlesUKUS.csv")

measles_UK <- measles_data %>% 
	filter(country=="UK") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases)) %>%
	filter(year < 1965)

measles_list <- measles_UK %>%
	split(as.character(.$loc))

nn <- names(tail(sort(sapply(measles_list, function(x) max(x$pop))), 10))

measles_list <- measles_list[nn] ## for teseting purposes

reconstruct_list <- measles_list %>%
	lapply(function(data) reconstruct_scam(data$cases, data$rec))

Zmat <- reconstruct_list %>%
	lapply("[[", "Z") %>%
	do.call(what="cbind")

rhomat <- reconstruct_list %>%
	lapply("[[", "rho") %>%
	do.call(what="cbind")

popmat <- measles_list %>%
	lapply("[[", "pop") %>%
	do.call(what="cbind") 

incmat <- measles_list %>%
	lapply("[[", "cases") %>%
	do.call(what="cbind")

## true (?) incidence
adjincmat <- round(incmat * rhomat)

distmat <- measles_list %>%
	lapply(filter, decimalYear==min(decimalYear)) %>%
	bind_rows(.id="city") %>%
	select(lon, lat) %>%
	dist(diag=TRUE, upper=TRUE) %>%
	as.matrix

invdist <- 1/distmat
diag(invdist) <- 0

## M[i, j] represents j to i movement
M <- popmat[1,] * invdist

x <- rep(1:52, 100)[1:(nrow(popmat)-1)]; k <- seq(0, 52, by=2)
BX <- cSplineDes(x,k)

standata <- list(
	ncity=ncol(popmat),
	nobs=nrow(popmat)-1,
	nbasis=ncol(BX),
	Inew=t(tail(adjincmat, -1)),
	Iprev=t(head(adjincmat, -1)),
	Zmat=t(head(Zmat, -1)),
	popmat=t(head(popmat, -1)),
	logpopmat=log(t(head(popmat, -1))),
	Bmat=t(BX),
	M=M/1e7
)
 
rt <- stanc(file="model.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

set.seed(101)
system.time(fit <- sampling(sm, data=standata, chains=1, iter=2000, thin=1))

save("nn", "fit", "standata", file="stan_UK.rda")
