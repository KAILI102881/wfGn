wFGN.eval <- function(H=NULL, n=1000, m=100, type="no", SNR=NULL, ndeps=c(1e-7,1e-2), noise=TRUE, pertype="per", minfun="qeta", weights=c(1,1,0), cluster=FALSE, plot=TRUE, sav=FALSE) {
	
	if(length(H)==0) {
		nloop <- 9; Hi <- 0; Hdelta <- 0.1
	} else {
		if((H>0)&(H<1)) {
			nloop <- 1; Hi <- H; Hdelta <- 0
		} else {
			stop("Invalid H")
		}
	}
	
	Hdata <- c(); Theta <- c(); Hstat <- c(); SNRdata <- c(); SNRstat <- c()
	for(iloop in (1:nloop)) {
		Hi <- Hi + Hdelta
		
		# replications of fGn
		theta1 <- c();	h <- c(); snr <- c()
		for(i in (1:m)) {	
			fGn <- perturbFGN(n, Hi, type=type, SNR=SNR, plot=FALSE)
			result <- wFGN(fGn, nloop=1, noise=noise, pertype=pertype, minfun=minfun,
					ndeps=ndeps, weights=weights, cluster=cluster, print.level=0)
			theta1[i] <- result$thetavector[1]
			h[i] <- result$thetavector[2]
			snr[i] <- result$SNR
		}
		
		# mean and sd of estimates
		Hmean <- mean(h)
		Hsd <- sd(h)
		MSE <- Hsd^2/n
		
		Hdata <- cbind(Hdata,h)
		Theta <- cbind(Theta,theta1)
		Hstat <- cbind(Hstat,list(H=Hi,Hmean=Hmean,Hsd=Hsd,MSE=MSE))
		
		SNRmean <- mean(snr)
		SNRsd <- sd(snr)
		SNR_MSE <- SNRsd^2/n
		
		SNRdata <- cbind(SNRdata,snr)
		SNRstat <- cbind(SNRstat,list(SNR=SNR,SNRmean=SNRmean,SNRsd=SNRsd,SNR_MSE=SNR_MSE))
	}
	
	if(plot) {
		Nchar <- as.character(n)
		Nleg <- paste(c("N= ", Nchar), collapse = " ")
		leg <- paste(c("Whittle estimator with ",Nleg))
		
		if(nloop==9) {
			boxplot(Hdata, main=leg,
					xlab="Hurst exponent",
					ylab="Estimated Hurst exponent", 
					names= seq(0.1,0.9,by=0.1))
		} else {
			boxplot(Hdata, main=leg,
					xlab="Hurst exponent",
					ylab="Estimated Hurst exponent", 
					names= c(H))
		}
		axis(2,seq(0.1,0.9,by=0.1))
	}
	drop(list(H=H,n=n,m=m,type=type,SNR=SNR,Hdata=Hdata,Theta=Theta,
					Hstat=Hstat,SNRdata=SNRdata,SNRstat=SNRstat))
	
	if(sav) {
		save(H,n,m,type,SNR,Hdata,Theta,Hstat,SNRdata,SNRstat,file = "wFGN.eval.RData")
	}
}