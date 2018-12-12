library(kernlab)
library(mgcv)
library(rstan)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("../R/reconstruct.R")
source("../R/summary.R")
source("../R/simulate.R")

load("stan_UK.rda")

measles_data <- read.csv("measlesUKUS.csv")

measles_UK <- measles_data %>% 
	filter(country=="UK") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases)) %>%
	filter(year < 1965)

measles_list <- measles_UK %>%
	split(as.character(.$loc))

reconstruct_list <- measles_list %>%
	lapply(function(data) reconstruct_gauss(data$cases, data$rec))

birthmat <- measles_list %>%
	lapply("[[", "rec") %>%
	do.call(what="cbind") %>%
	head(-1) %>%
	t

rhomat <- reconstruct_list %>%
	lapply("[[", "rho") %>%
	do.call(what="cbind") %>%
	t

decimalYear <- measles_list[[1]]$decimalYear

ext <- rstan::extract(fit)

x <- 1:52; k <- seq(0, 52, by=2)
BX <- cSplineDes(x,k)

psample <- seq(10, 250, by=10)

sumlist <- vector('list', length(psample))
caselist <- vector('list', length(psample))

set.seed(101)
for (j in 1:length(psample)) {
	i <- psample[j]
	betamat <- exp(ext$amat[i,,] %*% t(BX))
	
	alpha <- ext$alpha[i,]
	
	mixmat <- ext$m[i,,]
	
	popmat <- standata$popmat
	
	I0 <- standata$Iprev[,1]
	S0 <- round(ext$sbar[i,,1] + standata$Zmat[,1])
	
	sim <- simulate.sir(betamat, alpha, mixmat, I0, S0, popmat=popmat, birthmat=birthmat, rhomat=rhomat)
	
	sumlist[[j]] <- lapply(sim, function(x){
		cc <- x$C
		
		data.frame(
			city=1:nrow(cc),
			max=apply(cc, 1, max),
			zero=apply(cc, 1, function(y) sum(y==0)/ncol(cc)),
			period=apply(cc, 1, pfun)
		)
	}) %>%
		bind_rows(.id="sim")
	
	caselist[[j]] <- lapply(sim, function(x){
		cc <- x$C
		
		data.frame(
			city=rep(1:nrow(cc), each=ncol(cc)),
			tvec=rep(decimalYear, nrow(cc)),
			cases=c(t(cc))
		)
		
	}) %>%
		bind_rows(.id="sim")
}

summdf <- sumlist %>%
	bind_rows(.id="param") %>%
	mutate(city=factor(city, levels=1:40, labels=unique(measles_UK$loc))) %>%
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
			zero=sum(x==0)/nrow(x),
			period=pfun(x$cases)
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
	geom_smooth(aes(tvalue, mean), method="lm", se=FALSE, lty=2, fullrange=TRUE) +
	xlab("True maxima") +
	ylab("Simulated maxima")

g2 <- ggplot(combdf %>% filter(key=="zero")) +
	geom_point(aes(tvalue, mean)) +
	geom_errorbar(aes(tvalue, ymin=lwr, ymax=upr), width=0.0005) +
	geom_abline(intercept=0, slope=1) +
	geom_smooth(aes(tvalue, mean), method="lm", se=FALSE, lty=2, fullrange=TRUE) +
	xlab("True proportion of zeroes") +
	ylab("Simulated proportion of zeroes")

g3 <- ggplot(combdf %>% filter(key=="period")) +
	geom_point(aes(tvalue, mean)) +
	geom_errorbar(aes(tvalue, ymin=lwr, ymax=upr), width=0.0005) +
	geom_abline(intercept=0, slope=1) +
	geom_smooth(aes(tvalue, mean), method="lm", se=FALSE, lty=2, fullrange=TRUE) +
	xlab("True period") +
	ylab("Simulated period")

gtot <- arrangeGrob(g1, g2, g3, nrow=1)

ggsave("summary_UK.pdf", gtot, width=12, height=6)

casedf <- caselist %>%
	bind_rows(.id="param") %>%
	mutate(city=factor(city, levels=1:40, labels=unique(measles_UK$loc))) %>%
	group_by(city, tvec)

casesumm <- casedf %>%
	summarize(
		mean=mean(cases),
		lwr=quantile(cases, 0.025),
		upr=quantile(cases, 0.975)
	)

truecase <- measles_list %>%
	lapply(function(x){
		data.frame(
			city=unique(x$loc),
			cases=x$cases,
			tvec=decimalYear
		)
	}) %>%
	bind_rows %>%
	gather(key, value, -city, -tvec) %>%
	rename(tvalue=value)

gof <- merge(casesumm, truecase) %>%
	group_by(city) %>%
	summarize(
		gof=cor(mean, tvalue)^2
	) %>%
	arrange(-gof) %>%
	mutate(gof=formatC(gof, 3))

truecase$city <- factor(truecase$city, levels=gof$city)
casesumm$city <- factor(casesumm$city, levels=gof$city)

gcase <- ggplot(truecase) +
	geom_line(data=casesumm, aes(tvec, mean), col=2) +
	geom_ribbon(data=casesumm, aes(tvec, ymin=lwr, ymax=upr), fill=2, alpha=0.5) +
	geom_text(data=gof, aes(x=Inf, y=Inf, label=gof), hjust=1.05, vjust=1.2) +
	geom_point(aes(tvec, tvalue), size=0.1) +
	facet_wrap(~city, scale="free_y", nrow=8) +
	xlab("time") +
	ylab("reported cases") +
	theme(
		strip.background = element_blank(),
		panel.spacing = grid::unit(0, "cm")
	)

ggsave("simulate_UK.pdf", gcase, width=12, height=16)
save("caselist", "sumlist", file="simulate_UK.rda")
