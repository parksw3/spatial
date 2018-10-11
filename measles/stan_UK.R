library(mgcv)
library(dplyr)
library(scam)
library(rstan)
source("fitfun.R")

measles_data <- read.csv("measlesUKUS.csv")

measles_UK <- measles_data %>% 
	filter(country=="UK") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases))

measles_list <- measles_UK %>%
	split(as.character(.$loc))

measles_list <- measles_list[10:20] ## for teseting purposes

reconstruct_list <- measles_list %>%
	lapply(function(data) reconstruct(data$cases, data$rec))

Zmat <- reconstruct_list %>%
	lapply("[[", "Z") %>%
	do.call(what="cbind")

rhomat <- reconstruct_list %>%
	lapply("[[", "rho") %>%
	do.call(what="cbind")

popmat <- measles_list %>%
	lapply("[[", "pop") %>%
	do.call(what="cbind") %>%
	"["(1,) %>%
	rep(each=nrow(rhomat)) %>%
	matrix(nrow=nrow(rhomat))

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

x <- rep(1:26, 100)[1:(nrow(popmat)-1)]; k <- seq(0, 26, by=2) + 0.5
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
	M=M/1e8
)

rt <- stanc(file="model.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

set.seed(101) 
system.time(fit <- sampling(sm, data=standata, chains=1, iter=2000, thin=1))

save("fit", "standata", file="stan_UK.rda")

ext <- extract(fit)

ext$alpha[1000,]

plogis(ext$logit_sprop[1000,])

m <- ext$theta[1000] * standata$M
diag(m) <- 1 - colSums(m)

plot(c(standata$Inew[1,]), type="l")
lines(exp(ext$logIhat[1000,1,]), col=2)

x <- seq(0.5, 26.5, by=0.1)
BX2 <- cSplineDes(x,k)

plot(x, exp(c(apply(ext$mu_a, 2, median) %*% t(BX2))), type="l", ylim=c(20, 80))

for (i in 1:11) {
	lines(x, exp(c(ext$amat[1000,i,] %*% t(BX2))), type="l", col=2, lty=2)
}
