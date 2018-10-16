library(scam)
library(gamm4)
library(tsiR)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
source("fitfun.R")

load("analysis_UK.rda")

save <- FALSE

ff <- fitfun(coef(fit)[1], 0, 0, 0, distmat, adjincmat, Smat, popmat)

cmode <- lapply(colnames(popmat), function(x) {
	pp <- predict(ff$gamm4fit$gam,
				  newdata=data.frame(
				  	period=seq(1, 52, by=0.1),
				  	logI=0,
				  	off=0,
				  	region=x
				  ),
				  level=0.95,
				  se.fit=TRUE)
	
	data.frame(
		period=seq(1, 52, by=0.1),
		estimate=exp(pp$fit),
		lwr=exp(pp$fit - 1.96 * pp$se.fit),
		upr=exp(pp$fit + 1.96 * pp$se.fit),
		se=pp$se.fit,
		region=x
	)
}) %>%
	bind_rows

sumquad <- function(x) {sqrt(sum(c(x)^2))}

meanmode <- cmode %>%
	group_by(period) %>%
	summarise(estimate=mean(log(estimate)),
			  se=sumquad(se)) %>%
	mutate(
		upr=exp(estimate+1.96*se),
		lwr=exp(estimate-1.97*se),
		estimate=exp(estimate)
	)

alphamode <- data.frame(
	region=colnames(popmat),
	alpha=fixef(ff$gamm4fit$mer)[2] + ranef(ff$gamm4fit$mer)$region$logI
) %>%
	mutate(
		text=as.character(sprintf("%.3f", round(alpha, 3)))
	)

corder <- cmode %>%
	group_by(region) %>%
	summarize(
		cc=cor(estimate, meanmode$estimate)
	) %>%
	arrange(-cc)

cmode$region <- factor(cmode$region, level=corder$region)

gcond <- ggplot(cmode) + 
	geom_line(data=meanmode, aes(period, estimate)) +
	geom_line(aes(period, estimate), col="red", lty=2, lwd=1) + 
	geom_ribbon(aes(period, ymin=lwr, ymax=upr), alpha=0.1, fill="red", lty=2, col="red") +
	geom_text(data=alphamode, x=Inf, y=Inf, aes(label=text), hjust=1.05, vjust=1.25) +
	scale_x_continuous(expand=c(0,0)) +
	scale_y_continuous("transmission rate") +
	facet_wrap(~region, nrow=8) +
	theme(
		panel.grid = element_blank(),
		panel.spacing = grid::unit(0, "cm"),
		strip.background = element_blank()
	)

if (save) ggsave("conditional_mode_UK.pdf", gcond, width=8, height=12)
