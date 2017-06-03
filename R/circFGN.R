circFGN <- function(n, H, C=1, plot=FALSE) {
	
	# The below codes were copied from the package "longmemo".
	if(1 > (n <- as.integer(n))) stop("'n' must be a positive integer")
	H2 <- 2*H
	k <- 0:(n-1)
	autocov <- 0.5 * (abs(k-1)^H2 - 2*(k)^H2 + (k+1)^H2)
	
	m <- length(autocov)
	if(m < 2) stop("need at least 2 autocovariances")
	
	gk <- Re(fft(autocov[1:1 + c(0:(m - 1), if(m >= 3) (m - 2):1)],
					inverse = TRUE))
	
	if(any(gk < 0))
		stop("gk = Re(fft(autocov[1:n:2])) not all >= 0; invalid autocov[]")
	
	zr <- rnorm(m)
	zi <- rnorm(m-2)
	
	zr[c(1,m)] <- zr[c(1,m)]*sqrt(2)
	zr <- c(zr[1:m], zr[(m-1):2])
	zi <- c(0,zi,0,-zi[(m-2):1])
	
	z <- Re(fft((zr + 1i* zi) * sqrt(gk), inverse = TRUE))
	
	fGn <- ts(z[1:m] / (2*sqrt(m-1)))
	fGn <- fGn/sd(fGn)*sqrt(C)
	
	if(plot) {
		par(mfrow = c(1, 1))
		time <- (0:(n - 1))/n
		Nchar <- as.character(n)
		Nleg <- paste(c("N= ", Nchar), collapse = " ")
		Hchar <- as.character(round(H, 3))
		Hleg <- paste(c(", H=", Hchar), collapse = "")
		NHleg <- paste(c(Nleg, Hleg), collapse = "")
		leg <- paste(c(
						"Path of a fractional Gaussian noise ----- parameters",
						NHleg), collapse = " : ")
		plot(time, fGn, type = "l", main = leg)
	}	
	fGn
}