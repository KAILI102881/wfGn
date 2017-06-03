perturbFGN <- function(n, H, C=1, type="no", p=.005, SNR=NULL, plot=FALSE) {
	# Based on perturbFBM of the package dvfBm. Modified by W. You. Aug 2010.
	
	z <- circFGN(n, H, C=C, plot=FALSE)
	if(length(SNR)==0){
		zz <- z
	} else {
		switch(type,
				"no"={
					zz <- z
				},
				"WN"={ 
					s2B<-10^(-SNR/10)*C
					zz<-z+rnorm(n,mean=0,sqrt(s2B))
				},
				"AO"={
					ind<-which(rbinom(n-1,1,p)==1)
					if (length(ind)==0) ind<-sample(1:(n-1),1)
					nb<-length(ind)
					zz<-z
					s2Pert<- 10^(-SNR/10)*C
					pert<-sqrt(s2Pert)*rnorm(nb)
					zz[ind]<-zz[ind]+pert
				}
		)
	}
	
	if (plot) {
		txtMain<-paste("n=",n,", H=",H,sep="")
		if (type=="no") {
			par(mfrow=c(1,1))
			plot(zz,type="l",main=txtMain,xlab="Time",ylab="FGN")
		}
		if (type=="WN") {
			txtMain2<-paste("n=",n,", H=",H,", SNR=",SNR,sep="")
			par(mfrow=c(1,1))
			plot(zz,type="l",main=txtMain2,xlab="Time",ylab="noisy FGN")
		}
		if (type=="AO") {
			txtMain2<-paste("n=",n,", H=",H,", SNR=",SNR,sep="")
			par(mfrow=c(1,1))
			plot(zz,type="l",main=txtMain2,xlab="Time",ylab="FGN + outliers")
		}
		par(mfrow=c(1,1))
	}
	zz
}

