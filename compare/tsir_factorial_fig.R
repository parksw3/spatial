library(dplyr)

load("tsir_factorial.rda")

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

resdf <- reslist %>%
	bind_rows

bestgof <- which.min(resdf$ss2)
bestreg <- which.max(resdf$logLik)

pdf("tsir_factorial_fig.pdf", width=10, height=3.5)

par(mfrow=c(1, 3))

par(mar=c(5, 4, 4, 6) + 0.1)
plot(resdf$alpha, resdf$gof, type="l", ylim=c(0,0.85),
	 xlab="", ylab="", lty=2)
mtext("Goodness of fit",side=2,col="black",line=2.5)
mtext(expression(alpha),side=1,col="black",line=2.5)
lines(resdf$alpha, resdf$gof2, lty=1)
points(resdf$alpha[bestgof], resdf$gof2[bestgof], pch=16)
box()

par(new=TRUE)

plot(resdf$alpha, resdf$logLik, col="red", type="l", xlab="", ylab="", axes=FALSE, ylim=c(-210, -190))
mtext("Regression log-likelihood",side=4,col="red",line=4)
axis(4, col="red",col.axis="red",las=1)

points(resdf$alpha[bestreg], resdf$logLik[bestreg], pch=16, col=2)
legend(
	"topleft",
	legend=c("simulation", "simulation (init.fits)", "regression"),
	text.col=c("black", "black", "red"),
	lty=c(2, 1, 1),
	col=c(1, 1, 2)
)

title("a", adj=0, cex.main=2)

par(mar=c(5, 4, 4, 2) + 0.1)

plot(boston$time, fitlist[[bestgof]]$res$cases, 
	 xlab="", ylab="", pch=16, cex=0.5)
mtext("Biweekly incidence",side=2,col="black",line=2.5)
mtext("Time (years)",side=1,col="black",line=2.5)
lines(boston$time, fitlist2[[bestgof]]$res$V1)
lines(boston$time, fitlist2[[bestreg]]$res$V1, col=2)

title("b", adj=0, cex.main=2)


par(mar=c(5, 4, 4, 2) + 0.1)

plot(fitlist2[[bestgof]]$contact$beta * mean(boston$pop), type="l", ylim=c(10, 55),
	 xlab="",
	 ylab="")
mtext("Transmission rate",side=2,col="black",line=2.5)
mtext("Time (biweeks)",side=1,col="black",line=2.5)
polygon(c(1:26, 26:1), c(fitlist2[[bestgof]]$contact$betalow, rev(fitlist2[[bestgof]]$contact$betahigh)) * mean(boston$pop),
		col=rgb(0,0,0,0.3),
		border=NA)
lines(fitlist2[[bestreg]]$contact$beta * mean(boston$pop), col=2)
polygon(c(1:26, 26:1), c(fitlist2[[bestreg]]$contact$betalow, rev(fitlist2[[bestreg]]$contact$betahigh)) * mean(boston$pop),
		col=rgb(1,0,0,0.3),
		border=NA)

title("c", adj=0, cex.main=2)

dev.off()
