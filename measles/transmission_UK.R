library(mgcv)
library(rstan)
library(dplyr)
library(ggplot2); theme_set(theme_bw())

load("stan_UK.rda")

ext <- rstan::extract(fit)

x <- seq(0, 52, by=0.1); k <- seq(0, 52, by=2)
BX <- cSplineDes(x,k)

tfun <- function(coefmat) {
	tmat <- exp(apply(coefmat, 1, function(y) y %*% t(BX)))
	
	data.frame(
		period=x,
		mean=apply(tmat, 1, mean),
		lwr=apply(tmat, 1, quantile, 0.025),
		upr=apply(tmat, 1, quantile, 0.975)
	)
}

tcountry <- tfun(ext$mu_a)

tlist <- apply(ext$amat, 2, tfun)

names(tlist) <- rownames(standata$Inew)

tregion <- tlist %>%
	bind_rows(.id="city")

corder <- tregion %>%
	group_by(city) %>%
	summarize(cor=cor(mean, tcountry$mean)) %>%
	arrange(-cor) %>%
	select(city) %>%
	unlist %>%
	unname

tregion <- tregion %>%
	mutate(city=factor(city, level=corder))

alphadata <- apply(ext$alpha, 2, function(x){
	data.frame(t(sprintf("%.2f", c(
		mean=mean(x),
		lwr=quantile(x, 0.025),
		upr=quantile(x, 0.975)
	))))
}) %>%
	lapply(setNames, c("mean", "lwr", "upr")) %>%
	setNames(rownames(standata$Inew)) %>%
	bind_rows(.id="city") %>%
	mutate(
		text=paste0(mean, " (", lwr, ", ", upr, ")")
	) %>%
	mutate(city=factor(city, level=corder))

g1 <- ggplot(tregion) +
	geom_line(data=tcountry, aes(period, mean), lwd=1, col="grey") +
	geom_line(aes(period, mean), col=2, lty=1, lwd=1) +
	geom_ribbon(aes(period, ymin=lwr, ymax=upr), alpha=0.2, fill=2, col=2, lty=2) +
	geom_text(data=alphadata, x=Inf, y=Inf, aes(label=text), hjust=1.05, vjust=1.25) +
	facet_wrap(~city, nrow=5) +
	scale_x_continuous("biweek", expand=c(0,0)) +
	scale_y_continuous("transmission rate", limits=c(10, 80)) +
	theme(
		panel.grid = element_blank(),
		strip.background = element_blank(),
		panel.spacing = grid::unit(0, "cm")
	)

ggsave("transmission_UK.pdf", g1, width=16, height=16)
