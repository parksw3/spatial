library(kernlab)
library(tsiR)
library(dplyr)

rr <- read.csv("../measles/measlesUKUS.csv")

boston <- rr %>%
	filter(loc=="BOSTON") %>%
	filter(year >= 1920, year < 1940) %>%
	rename(
		time=decimalYear,
		cases=cases,
		pop=pop,
		births=rec
	)

alpha_vec <- seq(0.9, 1, by=0.0005)

basefit <- runtsir(boston,
				   Smean=seq(0.01*mean(boston$pop), 0.1*mean(boston$pop), length=50))	

fitlist <- fitlist2 <- reslist <- vector('list', length(alpha_vec))

for (i in 1:length(alpha_vec)) {
	print(i)
	
	alpha <- alpha_vec[i]
	
	bsfit <- runtsir(boston, regtype="user", userYhat=basefit$Yhat, alpha=alpha,
					 nsim=1,
					 Smean=seq(0.01*mean(boston$pop), 0.1*mean(boston$pop), length=50))
	
	bsfit2 <- runtsir(boston, regtype="user", userYhat=basefit$Yhat, alpha=alpha,
					  sbar=bsfit$sbar/mean(boston$pop),
					  nsim=1,
					  Smean=seq(0.01*mean(boston$pop), 0.1*mean(boston$pop), length=50),
					  inits.fit = TRUE,
					  nsample=50)
	
	fitlist[[i]] <- bsfit
	fitlist2[[i]] <- bsfit2
	
	reslist[[i]] <- data.frame(
		logLik=logLik(bsfit$glmfit)[[1]],
		gof=cor(boston$cases, bsfit$res[,1])^2,
		gof2=cor(boston$cases, bsfit2$res[,1])^2,
		ss=sum((boston$cases-bsfit$res[,1])^2),
		ss2=sum((boston$cases-bsfit2$res[,1])^2),
		alpha=alpha
	)
}

save("fitlist", "fitlist2", "reslist", file="tsir_factorial.rda")
