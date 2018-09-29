library(kernlab)
library(dplyr)
library(bbmle)
library(gamm4)
source("fitfun.R")

measles_data <- read.csv("measlesUKUS.csv")

measles_UK <- measles_data %>% 
	filter(country=="UK") %>%
	mutate(cases=ifelse(is.na(cases), 0, cases))

measles_list <- measles_UK %>%
	split(as.character(.$loc))

reconstruct_list <- measles_list %>%
	lapply(function(data) reconstruct(data$cases, data$rec))

Zmat <- reconstruct_list %>%
	lapply("[[", "Z") %>%
	do.call(what="cbind")

rhomat <- reconstruct_list %>%
	lapply("[[", "rho") %>%
	do.call(what="cbind")

popsize <- measles_list %>%
	lapply("[[", "pop") %>%
	do.call(what="cbind") %>%
	"["(1,)

## assumption
sbar <- 0.035 * popsize

Smat <- Zmat + matrix(rep(sbar, nrow(Zmat)), ncol=ncol(Zmat), byrow=TRUE)

any(Smat < 0)

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

fit <- mle2(nllfun,
	 start=list(log.theta=-20),
	 fixed=list(log.tau=0, log.rho=0, log.scale=0),
	 data=list(
	 	distmat=distmat,
	 	Imat=adjincmat,
	 	Smat=Smat,
	 	popsize=popsize,
	 	verbose=TRUE
	 )
)

fit <- fitfun(-20, 0, 0, 0, distmat, adjincmat, Smat, popsize)

cmode <- lapply(names(popsize), function(x) {
	pp <- predict(fit$gamm4fit$gam,
				  newdata=data.frame(
				  	period=seq(1, 26, by=0.1),
				  	logI=0,
				  	off=0,
				  	region=x,
				  	dummy=1
				  ),
				  level=0.95,
				  se.fit=TRUE)
	
	data.frame(
		period=seq(1, 26, by=0.1),
		estimate=exp(pp$fit),
		lwr=exp(pp$fit - 1.96 * pp$se.fit),
		upr=exp(pp$fit + 1.96 * pp$se.fit),
		region=x
	)
}) %>%
	bind_rows

ggplot(cmode) + 
	geom_line(aes(period, estimate)) + 
	geom_line(aes(period, upr), lty=2) + 
	geom_line(aes(period, lwr), lty=2) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous("transmission rate") +
	facet_wrap(~region, nrow=8) +
	theme(
		panel.grid = element_blank(),
		panel.spacing = grid::unit(0, "cm"),
		strip.background = element_blank()
	)
