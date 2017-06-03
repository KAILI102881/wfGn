wFGN <- function(x, istart=1, iend=length(x), nloop=1, init=c(0.55,0.01), ndeps=c(1e-7,1e-2), noise=TRUE, pertype="per", minfun="qeta", weights=c(1,1,0), cluster=FALSE, print.level=2) {
	# Author: Wonsang You 2010. This function is based on the S-PLUS codes 
	# originally developed by Jan Beran.
	
	require(multitaper)
	require(rgenoud)
	require(snow)
	require(Matching)
	
	# read data
	nmax <- length(x)
	n <- trunc((iend-istart+1)/nloop)
	nhalfm <- trunc((n-1)/2)
	
	# initialize h
	h <- init[1]
	eta <- c(h)
	
	# loop
	output <- c()
	thetavector <- c()
	i0 <- istart;
	for(iloop in (1:nloop)){
		h <- max(0.2, min(h,0.9)) # avoid extreme initial values
		eta[1] <- h
		i1 <- i0+n-1
		y <- x[i0:i1]
		
		# standardize data
		vary <- var(y)
		y <- (y-mean(y))/sqrt(var(y))
		
		# periodogram of data
		if(pertype=="taper"){
			yper <- spec.mtm(x,nw=4,k=8,plot=FALSE)$spec #length=padding+1
		}else{
			yper <- per(y)[2:(nhalfm+1)]
		}
		
		# find estimate
		if(noise) {
			InitVal<-init
			etaBnd<-c(0.001,0.999)
			snrBnd<-c(-10,10)
			DomainMat<-matrix(c(etaBnd,snrBnd), nrow = 2, ncol=2, byrow=TRUE)
			if(minfun=="lp") {
				result<-genoud(Lmin, nvars=2, starting.values=InitVal,
						max.generations=20, wait.generations=5, gradient.check=TRUE, 
						Domains=DomainMat, solution.tolerance=0.01, debug=TRUE,
						control=list(ndeps=ndeps), boundary.enforcement=2,
						cluster=cluster, print.level=print.level,
						n=n, yper=yper, pertype=pertype)		
			} else if(minfun=="csum"){
				result<-genoud(Cmin, nvars=2, starting.values=InitVal,
						max.generations=20, wait.generations=5, gradient.check=TRUE, 
						Domains=DomainMat, solution.tolerance=0.01, debug=TRUE,
						control=list(ndeps=ndeps), boundary.enforcement=2,
						cluster=cluster, print.level=print.level,
						n=n, yper=yper, pertype=pertype)
			} else if(minfun=="qeta"){
				result<-genoud(Qmin, nvars=2, starting.values=InitVal,
						max.generations=20, wait.generations=5, gradient.check=TRUE, 
						Domains=DomainMat, solution.tolerance=0.01, debug=TRUE,
						control=list(ndeps=ndeps), boundary.enforcement=2, 
						cluster=cluster, print.level=print.level,
						n=n, yper=yper, pertype=pertype)
			} else if(minfun=="combi"){
				result<-genoud(QLCmin, nvars=2, starting.values=InitVal,
						max.generations=20, wait.generations=5, gradient.check=TRUE, 
						Domains=DomainMat, solution.tolerance=0.01, debug=TRUE,
						control=list(ndeps=ndeps), boundary.enforcement=2, 
						cluster=cluster, print.level=print.level,
						n=n, yper=yper, pertype=pertype, weights=weights)
			}
			eta<-result$par[1]
			SNR<-result$par[2]
		} else {
			SNR <- NULL
			etatry <- eta
			if(minfun=="lp") {
				result <- optim(etatry, Lmin, NULL, method = "L-BFGS-B", 
						lower=1e-7, upper=1-1e-7, control=list(ndeps=ndeps[1]), 
						n=n, yper=yper, snr=SNR, pertype=pertype)
			} else if(minfun=="csum"){
				result <- optim(etatry, Cmin, NULL, method = "L-BFGS-B", 
						lower=1e-7, upper=1-1e-7, control=list(ndeps=ndeps[1]), 
						n=n, yper=yper, snr=SNR, pertype=pertype)
			} else if(minfun=="qeta"){
				result <- optim(etatry, Qmin, NULL, method = "L-BFGS-B", 
						lower=1e-7, upper=1-1e-7, control=list(ndeps=ndeps[1]), 
						n=n, yper=yper, snr=SNR, pertype=pertype)
			} else if(minfun=="combi"){
				result <- optim(etatry, QLCmin, NULL, method = "L-BFGS-B", 
						lower=1e-7, upper=1-1e-7, control=list(ndeps=ndeps[1]), 
						n=n, yper=yper, snr=SNR, pertype=pertype)
			}
			eta <- result$par
		}
		
		theta1 <- Qeta(eta,n,yper,SNR,pertype)$theta1
		theta <- c(theta1,eta)
		thetavector <- c(thetavector,theta)
		
		# calculate goodness of fit statistic
		Qresult <- Qeta(eta,n,yper,SNR,pertype)	
		
		# output
		M <- length(eta)
		SD <- CetaFGN(eta,snr=SNR)
		SD <- matrix(SD,ncol=M,nrow=M,byrow=TRUE)/n
		
		Hlow <- eta[1]-1.96*sqrt(SD[1,1])
		Hup  <- eta[1]+1.96*sqrt(SD[1,1])
		Hlow <- 0
		Hup <- 1
		fest <- Qresult$theta1*Qresult$fspec
		subseries <- list(thetavector=thetavector,
				Hlow=Hlow,Hup=Hup,SNR=SNR,fest=fest)
		output <- c(output, subseries)
		
		# next subseries
		i0 <- i0+n		
	}	
	drop(output)
}

per <- function(z) {	
	n <- length(z)
	(Mod(fft(z))^2/(2*pi*n)) [1:(n %/% 2 + 1)]
}

fspecFGN <- function(eta,m) {
	# parameters for the calculation of f
	h <- eta[1]
	nsum <- 200
	hh <- -2*h-1
	const <- 1/pi*sin(pi*h)*gamma(-hh)
	j <- 2*pi*c(0:nsum)
	
	# Fourier frequencies
	mhalfm <- trunc((m-1)/2)
	x <- 1:mhalfm
	x <- 2*pi/m*x
	
	# calculation of f at Fourier frequencies
	fspec <- matrix(0,mhalfm)
	for(i in seq(1:mhalfm)) {
		lambda <- x[i]
		fi <- abs(j+lambda)^hh+abs(j-lambda)^hh
		fi[1] <- fi[1]/2
		fi <- (1-cos(lambda))*const*fi
		fspec[i] <- sum(fi)
	}
	
	# adjusted spectrum (such that int(log(fspec))=0
	logfspec <- log(fspec)
	fint <- 2/(m)*sum(logfspec)
	theta1 <- exp(fint)
	fspec <- fspec/theta1
	drop(list(fspec=fspec,theta1=theta1))
}

fspecBFGN <- function(eta1, eta2, m) {
	
	# parameters for the calculation of f
	h1 <- eta1[1]
	h2 <- eta2[1]
	nsum <- 200
	hh <- -h1-h2-1
	const <- 1/pi*sin(pi/2*(h1+h2))*gamma(-hh)
	j <- 2*pi*c(0:nsum)
	
	# Fourier frequencies
	mhalfm <- trunc((m-1)/2)
	x <- 1:mhalfm
	x <- 2*pi/m*x
	
	# calculation of f at Fourier frequencies
	fspec <- matrix(0,mhalfm)
	for(i in seq(1:mhalfm)) {
		lambda <- x[i]
		fi <- abs(j+lambda)^hh+abs(j-lambda)^hh
		fi[1] <- fi[1]/2
		fi <- (1-cos(lambda))*const*fi
		fspec[i] <- sum(fi)
	}
	
	# adjusted spectrum such that int(log(fspec))=0
	logfspec <- log(fspec)
	fint <- 2/(m)*sum(logfspec)
	theta1 <- exp(fint)
	fspec <- fspec/theta1
	drop(list(fspec=fspec,theta1=theta1))
}

fspecPFGN <- function(eta, m, SNR=NULL) {
	
	if(length(SNR)>0){
		NvarBySvar <- 10^(-SNR/10)
	}else{
		NvarBySvar <- 0
	}
	
	Sresult1 <- fspecFGN(eta,m)
	fspec1 <- Sresult1$fspec
	
	Sresult2 <- fspecFGN(0.5,m)
	fspec2 <- Sresult2$fspec
	
	Sresult12 <- fspecBFGN(eta,0.5,m)
	fspec12 <- Sresult12$fspec
	
	if(length(SNR)>0){
		fspec <- fspec1*Sresult1$theta1+ fspec2*Sresult2$theta1*NvarBySvar+
				2*fspec12*Sresult12$theta1*sqrt(NvarBySvar)
	}else{
		fspec <- fspec1*Sresult1$theta1
	}
	
	logfspec <- log(fspec)
	fint <- 2/(m)*sum(logfspec)
	theta1 <- exp(fint)
	fspec <- fspec/theta1
	
	drop(list(fspec=fspec,theta1=theta1))
}

Qeta <- function(eta,n,yper,snr,pertype) {		
	h <- eta[1]
	
	if(pertype=="taper"){
		nn <- 2*(2^ceiling(log2(n))+1)+1
	}else{
		nn <- n
	}
	
	if(length(snr)==0) {
		fspec <- fspecFGN(eta,nn)
	} else {
		fspec  <- fspecPFGN(eta,nn,snr)
	}
	theta1 <- fspec$theta1
	fspec  <- fspec$fspec
	#logtheta1 <- log(theta1)
	
	yf <- yper/fspec
	yfyf <- yf^2
	A <- 2*(2*pi/nn)*sum(yfyf)
	B <- 2*(2*pi/nn)*sum(yf)
	Tn <- A/(B^2)
	z <- sqrt(nn)*(pi*Tn-1)/sqrt(2)
	pval <- 1-pnorm(z)
	theta1 <- B/(2*pi)
	fspec <- fspec

	#B <- logtheta1+B
	Qresult <- list(n=nn,h=h,
			eta=eta,A=A,B=B,Tn=Tn,z=z,pval=pval,
			theta1=theta1,fspec=fspec,value=B)
	drop(Qresult)
}

Qmin <- function(x,n,yper,pertype) {	
	# x[1]: etatry, x[2]: snr
	result <- Qeta(x[1],n,yper,x[2],pertype)$value
	drop(result)
}

Csum <- function(eta,n,yper,snr,pertype) {		
	h <- eta[1]
	
	if(pertype=="taper"){
		nn <- 2*(2^ceiling(log2(n))+1)+1
	}else{
		nn <- n
	}
	
	if(length(snr)==0) {
		fspec <- fspecFGN(eta,nn)
	} else {
		fspec  <- fspecPFGN(eta,nn,snr)
	}
	theta1 <- fspec$theta1
	fspec  <- theta1*fspec$fspec
	
	mhalfm <- trunc((nn-1)/2)
	cyper <- cumsum(yper)/sum(yper)
	cfspec <- cumsum(fspec)/sum(fspec)
	sqdiff <- (cyper-cfspec)^2
	sqdiffm <- sum(sqdiff)/mhalfm
	
	yf <- yper/fspec
	yfyf <- yf^2
	A <- 2*(2*pi/nn)*sum(yfyf)
	B <- 2*(2*pi/nn)*sum(yf)
	Tn <- A/(B^2)
	z <- sqrt(nn)*(pi*Tn-1)/sqrt(2)
	pval <- 1-pnorm(z)
	theta1 <- B/(2*pi)
	fspec <- fspec
	
	Qresult <- list(n=nn,h=h,
			eta=eta,A=A,B=B,Tn=Tn,z=z,pval=pval,
			theta1=theta1,fspec=fspec,value=sqdiffm)
	drop(Qresult)
}

Cmin <- function(x,n,yper,pertype) {	
	# x[1]: etatry, x[2]: snr
	result <- Csum(x[1],n,yper,x[2],pertype)$value
	drop(result)
}

Lpvar <- function(eta,n,yper,snr,pertype) {		
	h <- eta[1]
	
	if(pertype=="taper"){
		nn <- 2*(2^ceiling(log2(n))+1)+1
	}else{
		nn <- n
	}
	
	if(length(snr)==0) {
		fspec <- fspecFGN(eta,nn)
	} else {
		fspec  <- fspecPFGN(eta,nn,snr)
	}
	theta1 <- fspec$theta1
	fspec  <- fspec$fspec
	
	epsilon<-log(yper/(fspec*theta1))+0.57721
	
	yf <- yper/fspec
	yfyf <- yf^2
	A <- 2*(2*pi/nn)*sum(yfyf)
	B <- 2*(2*pi/nn)*sum(yf)
	Tn <- A/(B^2)
	z <- sqrt(nn)*(pi*Tn-1)/sqrt(2)
	pval <- 1-pnorm(z)
	theta1 <- B/(2*pi)
	fspec <- fspec
	
	Qresult <- list(n=nn,h=h,
			eta=eta,A=A,B=B,Tn=Tn,z=z,pval=pval,
			theta1=theta1,fspec=fspec, value=var(epsilon))
	drop(Qresult)
}

Lmin <- function(x,n,yper,pertype) {	
	# x[1]: etatry, x[2]: snr
	result <- Lpvar(x[1],n,yper,x[2],pertype)$value
	drop(result)
}

QLCfun <- function(eta,n,yper,snr,pertype,weights) {		
	h <- eta[1]
	
	if(pertype=="taper"){
		nn <- 2*(2^ceiling(log2(n))+1)+1
	}else{
		nn <- n
	}
	
	if(length(snr)==0) {
		fspec <- fspecFGN(eta,nn)
	} else {
		fspec  <- fspecPFGN(eta,nn,snr)
	}
	theta1 <- fspec$theta1
	fspec  <- fspec$fspec
	
	if(weights[2]>0){
		epsilon <- log(yper/(fspec*theta1))+0.57721
		vareps <- var(epsilon)
	}else{
		vareps <- 0
	}
	if(weights[3]>0){
		mhalfm <- trunc((nn-1)/2)
		cyper <- cumsum(yper)/sum(yper)
		cfspec <- cumsum(fspec)/sum(fspec)
		sqdiff <- (cyper-cfspec)^2
		sqdiffm <- sum(sqdiff)/mhalfm
	}else{
		sqdiffm <- 0
	}
	
	yf <- yper/fspec
	yfyf <- yf^2
	A <- 2*(2*pi/nn)*sum(yfyf)
	B <- 2*(2*pi/nn)*sum(yf)
	Tn <- A/(B^2)
	z <- sqrt(nn)*(pi*Tn-1)/sqrt(2)
	pval <- 1-pnorm(z)
	theta1 <- B/(2*pi)
	fspec <- fspec
	
	value=weights[1]*B+weights[2]*vareps+weights[3]*sqdiffm
	
	Qresult <- list(n=nn,h=h,
			eta=eta,A=A,B=B,Tn=Tn,z=z,pval=pval,
			theta1=theta1,fspec=fspec, value=value)
	drop(Qresult)
}

QLCmin <- function(x,n,yper,pertype,weights) {	
	# x[1]: etatry, x[2]: snr
	result <- QLCfun(x[1],n,yper,x[2],pertype,weights)$value
	drop(result)
}

CetaFGN <- function(eta,snr=NULL) {		
	M <- length(eta)
	
	# size of steps in Riemann sum: 2*pi/m
	m <- 10000
	mhalfm <- trunc((m-1)/2)
	
	# size of delta for numerical calculation of derivative
	delta <- 1e-9
	
	# partial derivatives of log f (at each Fourier frequency)
	lf <- matrix(rep(1,M*mhalfm),ncol=M,nrow=mhalfm)
	if(length(snr)==0) {
		f0 <- fspecFGN(eta,m)$fspec
		for(j in (1:M)) {
			etaj <- eta
			etaj[j] <- etaj[j]+delta
			fj <- fspecFGN(etaj,m)$fspec
			lf[,j] <- log(fj/f0)/delta
		}	
	} else {
		f0 <- fspecPFGN(eta,m,snr)$fspec
		for(j in (1:M)) {
			etaj <- eta
			etaj[j] <- etaj[j]+delta
			fj <- fspecPFGN(etaj,m,snr)$fspec
			lf[,j] <- log(fj/f0)/delta
		}	
	}
	
	# Calculate D
	Djl <- matrix(rep(1,M*M),ncol=M,nrow=M)
	for(j in (1:M)) {
		for(l in (1:M)) {
			Djl[j,l] <- 2*2*pi/m*sum(lf[,j]*lf[,l])
		}
	}
	drop(matrix(4*pi*solve(Djl),ncol=M,nrow=M,byrow=TRUE))
}
