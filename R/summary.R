pfun <- function(x, span=5) {
	ss <- spectrum(x, log="no", span=span, plot=FALSE)
	
	1/ss$freq[which.max(ss$spec)]/26
}