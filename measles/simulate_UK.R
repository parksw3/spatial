library(rstan)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("fitfun.R")
source("simulate.R")

load("stan_UK.rda")

measles_data <- read.csv("measlesUKUS.csv")

measles_UK <- measles_data %>% 
	filter(country=="UK") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases))

measles_list <- measles_UK %>%
	split(as.character(.$loc))

nn <- names(tail(sort(sapply(measles_list, function(x) max(x$pop))), 10))

measles_list <- measles_list[nn] ## for teseting purposes

reconstruct_list <- measles_list %>%
	lapply(function(data) reconstruct(data$cases, data$rec))

birthmat <- measles_list %>%
	lapply("[[", "rec") %>%
	do.call(what="cbind") %>%
	head(-1) %>%
	t

rhomat <- reconstruct_list %>%
	lapply("[[", "rho") %>%
	do.call(what="cbind") %>%
	t

ext <- rstan::extract(fit)

x <- 1:52; k <- seq(0, 52, by=2) + 0.5
BX <- cSplineDes(x,k)

psample <- seq(10, 1000, by=10)

sumlist <- vector('list', 100)

for (j in 1:100) {
	i <- psample[j]
	tmat <- exp(ext$amat[i,,] %*% t(BX))
	
	alpha <- ext$alpha[i,]
	
	mixmat <- ext$m[i,,]
	
	popmat <- standata$popmat
	
	I0 <- standata$Iprev[,1]
	S0 <- round(ext$sbar[i,,1] + standata$Zmat[,1])
	
	set.seed(101)
	sim <- simfun(tmat, alpha, mixmat, I0, S0, popmat=popmat, birthmat=birthmat, rhomat=rhomat)
	
	sumlist[[j]] <- lapply(sim, function(x){
		cc <- x$C
		
		data.frame(
			city=nn,
			max=apply(cc, 1, max),
			zero=apply(cc, 1, function(y) sum(y==0)/ncol(cc))
		)
	}) %>%
		bind_rows(.id="sim")
}

summdf <- sumlist %>%
	bind_rows(.id="param") %>%
	gather(key, value, -city, -sim, -param) %>%
	group_by(key, city) %>%
	summarize(
		mean=mean(value),
		lwr=quantile(value, 0.025),
		upr=quantile(value, 0.975)
	)

truesumm <- measles_list %>%
	lapply(function(x){
		data.frame(
			city=unique(x$loc),
			max=max(x$cases),
			zero=sum(x==0)/nrow(x)
		)
	}) %>%
	bind_rows(.id="city") %>%
	gather(key, value, -city) %>%
	rename(tvalue=value)

combdf <- merge(summdf, truesumm) %>%
	mutate(tvalue=as.numeric(as.character(tvalue)))

g1 <- ggplot(combdf %>% filter(key=="max")) +
	geom_point(aes(tvalue, mean)) +
	geom_errorbar(aes(tvalue, ymin=lwr, ymax=upr), width=0.1) +
	geom_abline(intercept=0, slope=1) +
	xlab("True maxima") +
	ylab("Simulated maxima")

g2 <- ggplot(combdf %>% filter(key=="zero")) +
	geom_point(aes(tvalue, mean)) +
	geom_errorbar(aes(tvalue, ymin=lwr, ymax=upr), width=0.0005) +
	geom_abline(intercept=0, slope=1) +
	xlab("True proportion of zeroes") +
	ylab("Simulated proportion of zeroes")

gtot <- arrangeGrob(g1, g2, nrow=1)

ggsave("simulate_UK.pdf", gtot, width=8, height=6)
