library(scam)
library(dplyr)
library(bbmle)
library(gamm4)
library(tsiR)
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

popmat <- measles_list %>%
	lapply("[[", "pop") %>%
	do.call(what="cbind")

birthmat <- measles_list %>%
	lapply("[[", "rec") %>%
	do.call(what="cbind")

## assumption
sbar <- 0.035 * popmat

Smat <- Zmat + sbar

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
	 	popmat=popmat,
	 	verbose=TRUE
	 )
)

save("Zmat", "popmat", "birthmat", "Smat", "rhomat", "incmat", "adjincmat", 
	 "distmat", "fit", file="analysis_UK.rda")
