Umbrella2d <- function(datalist,vol=0.85,NTC=NULL,Plate=T,Rotate=F){
	# datalist must be a list of dataframes.
	# 	Each dataframe (or matrix) should have as first 2 columns fluorescence intensities of 2 channels
	# 	an optional column with prior cluster identification named "Cluster" may be added
	# vol is the volume of one droplet / partition in nL (0.89 for Bio-Rad's QX-100)
	# NTC is a vector of names indicating which elements from datalist are NTC's
	# 	The names in NTC must match the names of items in the datalist.
	# Plate: if true, all plate information may be used to calculate cluster centres (default)
	# 	if false, only cluster information (if given) and well data are used.
	# Rotate: if false, no data rotation is performed. (default)
	#	if true, data gets rotated based on a default, well-based algorithm.
	#	if "Plate", data gets rotated based on data from complete plate
	# 	if a name of a well is given, all wells are rotated based on data from this well.

require(MASS)
require(mgcv)
require(robust)
require(modeest)
require(OrdMonReg)

########################
######   Proc 1   ######
########################

### function to find the cluster centres
	# search center of double null and single positive clusters

modeRetFin=function(data){

# method based on the initial clusters (default)
data1 <- data[data$Cluster==nullclus[1],]
data2 <- data[data$Cluster==nullclus[2],]
data3 <- data[data$Cluster==nullclus[3],]
	# data in each cluster
clusn <- c(nrow(data1),nrow(data2),nrow(data3))/nrow(data)
	# percentage in each cluster
hlp <- rep(NA,6)
if(clusn[1]>0){hlp[1:2] <- c(mlv(data1[,1],method="asselin")$M,mlv(data1[,2],method="asselin")$M)}
if(clusn[2]>0){hlp[3:4] <- c(mlv(data2[,1],method="asselin")$M,mlv(data2[,2],method="asselin")$M)}
if(clusn[3]>0){hlp[5:6] <- c(mlv(data3[,1],method="asselin")$M,mlv(data3[,2],method="asselin")$M)}
names(hlp)=c("xMode00","yMode00","xMode10","yMode10","xMode01","yMode01")
	# identify each of the modes

# method based on single well to check for problems
cutx=mean(quantile(data[,1],c(0.01,0.999)))
cuty=mean(quantile(data[,2],c(0.01,0.999)))
	# arbitrary cutpoint
data1w <- data[data[,1]<cutx & data[,2]<cuty,]
data2w <- data[data[,1]>=cutx & data[,2]<cuty,]
data3w <- data[data[,1]<cutx & data[,2]>=cuty,]
	# data in each cluster
clusnw <- c(nrow(data1w),nrow(data2w),nrow(data3w))/nrow(data)
	# percentage in each cluster
hlp2 <- rep(NA,6)
if(clusnw[1]>0){hlp2[1:2] <- c(mlv(data1w[,1],method="asselin")$M,mlv(data1w[,2],method="asselin")$M)}
if(clusnw[2]>0){hlp2[3:4] <- c(mlv(data2w[,1],method="asselin")$M,mlv(data2w[,2],method="asselin")$M)}
if(clusnw[3]>0){hlp2[5:6] <- c(mlv(data3w[,1],method="asselin")$M,mlv(data3w[,2],method="asselin")$M)}
names(hlp2)=c("xMode00","yMode00","xMode10","yMode10","xMode01","yMode01")

hlp[is.na(hlp)] <- hlp2[is.na(hlp)]
hlp2[is.na(hlp2)] <- hlp[is.na(hlp2)]
	# replace NA's by the estimated version from the other method

moddiff <- hlp-hlp2
clusdiff <- rep(c(hlp[3]-hlp[1],hlp[6]-hlp2[2]),3)
suppressWarnings(diffratio <- moddiff/clusdiff)

diffc <- F
if(max(abs(diffratio),na.rm=T)>0.1){diffc <- T}
	# T/F whether the difference between estimated centres is large

if(Plate==F & clinm!="input"){hlp <- hlp2
					clusn <- clusnw}
	# No plate information used!

return(list(center=hlp,clusp=clusn,diffc=diffc))
}

########################
######   Proc 2   ######
########################

### function to check centres based on plate information

checkmod1 <- function(modsds){
	modout <- modsds
		# give warning when a cluster (centre) is not defined.
		# don't warn if it's the single positive of an NTC.
	if(is.na(modsds$center[1])|is.na(modsds$center[2])){
		modout$center[1:2] <- Cmat[1:2,1]
		cat("      Warning: (0,0) cluster in",modsds$name,"not defined.\n")}
	if(is.na(modsds$center[3])|is.na(modsds$center[4])){
		modout$center[3:4] <- modout$center[1:2]+Cmat[3:4,1]
		suppressWarnings(if(max(NTC == modsds$name)!=1){cat("      Warning: (1,0) cluster in",modsds$name,"not defined.\n")})}
	if(is.na(modsds$center[5])|is.na(modsds$center[6])){
		modout$center[5:6] <- modout$center[1:2]+Cmat[5:6,1]
		suppressWarnings(if(max(NTC == modsds$name)!=1){cat("      Warning: (0,1) cluster in",modsds$name,"not defined.\n")})}
	modd <- modout$center-c(0,0,rep(modout$center[1:2],2))
		# centres/distances for the specific well
	dists <- abs(modd-Cmat[,1])/Cmat[,2]
		# relative difference from plate median (number of robust sd's)
	dists[is.na(dists)] <- 0
	eps <- .Machine$double.eps
	if((dists[1]>5) | (dists[2]>5))
		{cat("      Warning: center (0,0) cluster in",modsds$name,"may deviate from other centers in plate. Please check input.\n")}
		# never replace double negative, only give warning.
	if( ((dists[3]+dists[4])/(modout$clusp[2]+eps) >1000) & (modout$clusp[2] < 0.01) )
		{modout$center[3] <- modout$center[1]+Cmat[3,1]
		modout$center[4] <- modout$center[2]+Cmat[4,1]
		suppressWarnings(if(max(NTC == modsds$name)!=1){cat("      Warning: center (1,0) cluster in",modsds$name,"may deviate differently from its (0,0) cluster center compared to other wells of the plate.  Plate information used instead.\n")})}
	if( ((dists[5]+dists[6])/(modout$clusp[3]+eps) >1000) & (modout$clusp[3] < 0.01) )
		{modout$center[5] <- modout$center[1]+Cmat[5,1]
		modout$center[6] <- modout$center[2]+Cmat[6,1]
		suppressWarnings(if(max(NTC == modsds$name)!=1){cat("      Warning: center (0,1) cluster in",modsds$name,"may deviate differently from its (0,0) cluster center compared to other wells of the plate.  Plate information used instead.\n")})}
		# replace single positive if limited data and difference large
	return(modout)
		}

########################
######   Proc 3   ######
########################

### function to rotate data

Rotfunc <- function(data2d,mods){
	x00 <- mods$center[1]
	y00 <- mods$center[2]
	x10 <- mods$center[3]
	y10 <- mods$center[4]
	x01 <- mods$center[5]
	y01 <- mods$center[6]
		# cluster centres
	ricox <- (y10-y00)/(x10-x00)
	ricoy <- (x01-x00)/(y01-y00)
		# slopes of connecting lines
	if(is.na(ricox)){ricox <- 0}
	if(is.na(ricoy)){ricoy <- 0}
		# 0 means a channel remains the same after rotation
	denom <- 1-ricox*ricoy
	datax <- (data2d[,1]-ricoy*data2d[,2]+y00*ricoy-x00)/denom+mods$center[1]
	datay <- (data2d[,2]-ricox*data2d[,1]+x00*ricox-y00)/denom+mods$center[2]
		# rotated data + original double negative center added
return(data.frame(datax,datay,data2d[,3]))
}

########################
######   Proc 4   ######
########################

# Fit a spline on the NTC

NTCfunc <- function(data){
	data <- sort(data)
	sigma <- huber(data)$s
		# robust deviation
	data <- data[(data > -10*sigma) & (data < 10*sigma)]
		# use only data within 10 deviations
	ntcbrks <- seq(-11*sigma,11*sigma,sigma/10)
		# breaks
	ntcfreq <- hist(data,breaks=ntcbrks,plot = FALSE)$counts
		# create table with frequencies (include zero's)
	corfac <- length(data)*sigma/10
		# correction to transform absolute number into density
return(list(NTC=ntcfreq,NTCden=corfac,NTCsig=sigma))
}


########################
######   Proc 5   ######
########################

# Fit spline and normal null distribution (2nd order polynomial) on the data
# back-up procedure for no NTC's

Locfunc1 <- function(datafit,mu,datapred=NULL){
	# datafit should be a datavector
	# mu should be the null mode
	# datapred is the prediction dataset, if different from the input data

if(is.null(datapred)){datapred <- datafit}
fdr <- NULL
nullfit <- NULL
datasplf <- NULL
mufit <- NULL
sigfit <-NULL
pifit <- NULL

if(length(datafit)>3){
### estimate sigma

Xnull <- datafit[datafit<mean(quantile(datafit,c(0.001,0.999)))]
sigma <- huber(Xnull)$s


### spline fit

sdss <- (max(datafit)-min(datafit))/sigma+2

spbins <- ceiling(sdss*10)
spbin <- seq(min(datafit)-sigma,max(datafit)+sigma,length=spbins+1)
spmids <- (spbin[-1]+spbin[-length(spbin)])/2
	# bin intervals
databin <- findInterval(datafit,spbin)
spfreq <- table(c(databin,1:spbins))-rep(1,spbins)
fmod <- gam(spfreq ~ s(spmids,k=min(ceiling(sdss),100)), family=poisson)
f <- predict(fmod,newdata=data.frame(spmids=datapred),type="response")

spbinw <- (max(datafit)+2*sigma-min(datafit))/spbins
datasplf <- f/length(datafit)/spbinw

sppred <- spmids[(spmids < mu+sigma) & (spmids > mu-2*sigma)]
spresp <- predict(fmod,newdata=data.frame(spmids=sppred),type="response")

sppredl <- sppred-mu
sppredq <- sppredl^2
sprespl <- log(spresp)

try(spqmodp <- lm(sprespl ~ sppredl+sppredq))
try(sigfit <- sqrt(-1/2/spqmodp$coef[3]))
try(mufit <- -spqmodp$coef[2]/(2*spqmodp$coef[3]))
try(pifit <- exp(spqmodp$coef[1]-log(length(datafit))+(mufit^2/sigfit^2+log(2*pi*sigfit^2))/2)/spbinw)
}

else{cat("      Warning: not enough data to fit null distribution\n")}

try(nullfit <- dnorm(datapred,mufit+mu,sigfit)*pifit)
try(nullfit[datapred < mufit+mu-sigfit] <- datasplf[datapred < mufit+mu-sigfit])
res <- c(mufit+mu,sigfit,pifit)
try(fdr <- nullfit/datasplf)
return(list(fdr=fdr,datanull=nullfit,datafit=datasplf,pifits=pifit,res=res))
}



########################
######   Proc 6   ######
########################

# Fit NTC on the data

Locfunc2 <- function(datafit,mu,NTC,datapred=NULL){
	# data should be a datavector
	# mu should be the null mode
	# NTC should be a list with NTC model:
	# NTC$NTC are the NTC frequencies
	# NTC$NTCden is the correction to get a density
	# NTC$NTCsig is the robust deviation (mad) of the NTC
	
	if(is.null(datapred)){datapred <- datafit}
	fdr <- NULL
	fmarg <- NULL
	f0 <- NULL
	pi0 <- NULL
	prob <- NULL
	sigma <- NTC$NTCsig
### spline fit

try(if(length(datafit)>3){
	data <- datafit-mu
	datap <- datapred-mu
	datmin <- min(floor(min(data)/sigma),-11)
	datmax <- max(ceiling(max(data)/sigma),11)
	datbrks <- seq(datmin*sigma,datmax*sigma,sigma/10)
		# breaks
	datmids <- (datbrks[-1]+datbrks[-length(datbrks)])/2
		# bin mids
	datfreq <- hist(data,breaks=datbrks,plot = FALSE)$counts
		# create table with frequencies (include zero's)
	ntcfreq <- c(rep(0,(-11-datmin)*10),NTC$NTC,rep(0,(datmax-11)*10))
	NTCden <- NTC$NTCden
	datden <- length(data)*sigma/10
		# corrections to transform counts into densities

	dataSam <- rbind(data.frame(mids=datmids,counts=ntcfreq,label="ntc",den=NTCden),data.frame(mids=datmids, counts=datfreq,label="well",den=datden))
	dataSam$dummy <- as.double(dataSam$label)-1
	eps <- -ceiling(log10(sigma/100))
	if(sum(ntcfreq>0)>3){
		midHlp.ntc <- datmids[ntcfreq>0]
		knots.ntc <- quantile(midHlp.ntc,seq(0,1,1/(min(10,length(midHlp.ntc)-3))))
			# 11 NTC knots (or equal to number of filled bins -2 if lower)
		}else{knots.ntc <- NULL}
	if(sum((datfreq>0)&(datmids>2*sigma)&(datmids<11*sigma))>1){
		midHlp.tail <- datmids[(datfreq>0)&(datmids>2*sigma)&(datmids<11*sigma)]
		knots.tail <- quantile(midHlp.tail,seq(0,1,1/(min(5,length(midHlp.tail)-1))))
			# 6 tail knots (or equal to number of filled bins if lower)
		}else{knots.tail <- NULL}
	if(sum((datfreq>0)&(datmids>11*sigma))>1){
		midHlp.pos <- datmids[(datfreq>0)&(datmids>11*sigma)]
		knots.pos <- quantile(midHlp.pos,seq(0,1,1/(min(5,length(midHlp.pos)-1))))
			# 6 pos knots (or equal to number of filled bins if lower)
		}else{knots.pos <- NULL}
	knots <- sort(unique(round(c(knots.ntc,knots.tail,knots.pos,0,min(datmids[datfreq>0])),eps)))
		# 1 knot at 0
		# 1 knot at minimum real data
		# 25 knots, possibly not unique, only leave unique knots
		# round to make sure knots do not differ by very small values.
	nknot <- length(knots)
	fitHlp <- gam(counts~s(mids,k=nknot)+s(mids,by=dummy,k=nknot)+dummy,knots=list(mids=knots),family=poisson,data=dataSam,offset=I(log(den)))

	dataPred0 <- data.frame(mids=rep(seq(-3,2,0.1),2)*sigma,dummy=rep(c(0,1),each=51))
	pi0=exp(diff(predict(fitHlp,dataPred0),51))
	prob <- mean(pi0[21:31])

	fmarg <- predict(fitHlp,data.frame(mids=datap,dummy=1),type="response")
	f0 <- prob*predict(fitHlp,data.frame(mids=datap,dummy=0),type="response")
	fdr <- f0/fmarg
	fdr[datap < -sigma] <- 1
}
else{cat("      Warning: not enough data to fit null distribution\n")})

res <- c(mu,NTC$NTCsig,prob)
return(list(fdr=fdr,datafit=fmarg,datanull=f0,pifits=pi0,res=res))
}


########################
######   Proc 7   ######
########################

# Shell to run proc 5, 6

Shellfunc <- function(data2d,modout){
	time1 <- Sys.time()
	Locoutx <- list(fdr=NULL,datanull=NULL,datafit=NULL,pifits=NULL,res=NULL)
	Locouty <- list(fdr=NULL,datanull=NULL,datafit=NULL,pifits=NULL,res=NULL)
	Locoutx0 <- list(fdr=NULL,datanull=NULL,datafit=NULL,pifits=NULL,res=NULL)
	Locouty0 <- list(fdr=NULL,datanull=NULL,datafit=NULL,pifits=NULL,res=NULL)
	Locoutx1 <- list(fdr=NULL,datanull=NULL,datafit=NULL,pifits=NULL,res=NULL)
	Locouty1 <- list(fdr=NULL,datanull=NULL,datafit=NULL,pifits=NULL,res=NULL)
		# create empty objects in case something goes wrong
	datax <- data2d[,1]
	datay <- data2d[,2]
	rax <- rank(data2d[,1],ties.method="first")
	ray <- rank(data2d[,2],ties.method="first")
	orx <- order(data2d[,1])
	ory <- order(data2d[,2])
	dataxs <- sort(datax)
	datays <- sort(datay)
		# univariate channel data and sort objects
	probx <- rep(NA,length(datax))
	proby <- rep(NA,length(datay))
	problistx <- list()
	problisty <- list()
	probrx <- probx
	probry <- proby
	probciLx <- probx
	probciUx <- probx
	probciLy <- proby
	probciUy <- proby
	sigdifx <- NA
	sigdify <- NA
		# create empty objects to store results

Ttime1 <- Sys.time()

# initial fit on full data
	if((!is.null(NTCxN))&(!is.null(NTCyN))){
		try(Locoutx <- lapply(as.list(1:length(NTCxN)),function(i){Locfunc2(dataxs,modout$center[1],NTCxN[[i]])}))
		try(Locouty <- lapply(as.list(1:length(NTCyN)),function(i){Locfunc2(datays,modout$center[2],NTCyN[[i]])}))
			# if NTC's are present, use Locfunc2 for fitting
			# fit on sorted data
		try(problistx <- lapply(Locoutx,function(x){x$fdr}))
		try(problisty <- lapply(Locouty,function(x){x$fdr}))
			# store posterior probabilities
	Ttime2 <- Sys.time()
		try(probrx <- matrix(unlist(problistx),ncol=length(NTCxN)))
		try(probry <- matrix(unlist(problisty),ncol=length(NTCyN)))
			# create matrix of probabilities with a column per NTC
		try(probmx <- rowMeans(probrx))
		try(probmy <- rowMeans(probry))
			# take average probability over the NTC's
		try(probx <- BoundedAntiMean(probmx,w=rep(1,length(datax)),a=rep(0,length(datax)),b=rep(1,length(datax)))[rax])
		try(proby <- BoundedAntiMean(probmy,w=rep(1,length(datay)),a=rep(0,length(datay)),b=rep(1,length(datay)))[ray])
			# bounded antitonic regression on the probabilities
			# sorted results are reordered using rax and ray (original ranks)
		} # end NTC procedure
	else{ try(Locoutx <- Locfunc1(dataxs,modout$center[1]))
		try(Locouty <- Locfunc1(datays,modout$center[2]))
			# if no NTC's are present, use Locfunc1 for fitting
		try(probrx <- Locoutx$fdr)
		try(probry <- Locouty$fdr)
			# store posterior probabilities
	Ttime2 <- Sys.time()
		try(probx <- BoundedAntiMean(probrx,rep(1,length(datax)),a=rep(0,length(datax)),b=rep(1,length(datax)))[rax])
		try(proby <- BoundedAntiMean(probry,rep(1,length(datay)),a=rep(0,length(datay)),b=rep(1,length(datay)))[ray])
			# bounded antitonic regression on the probabilities
			# sorted results are reordered using rax & ray (original ranks)
		mux <- Locoutx$res[1]
		sigx <- Locoutx$res[2]
		pix <- Locoutx$res[3]
		muy <- Locouty$res[1]
		sigy <- Locouty$res[2]
		piy <- Locouty$res[3]
		} # end back-up non-NTC procedure

		# create initial estimates and objects

Ttime3 <- Sys.time()
Ttime4 <- list()
for(i in 1:4){
Ttime4[[i]] <- Sys.time()
}

# update fit specific for negative/positive on other channel
# NTC version fit
	if((!is.null(NTCxN))&(!is.null(NTCyN))){
		try(datax0 <- datax[proby >= 0.05])
		try(datax1 <- datax[proby < 0.05])
		try(datay0 <- datay[probx >= 0.05])
		try(datay1 <- datay[probx < 0.05])
			# data split over other channel
		try(rax0 <- rank(datax0,ties.method="first"))
		try(rax1 <- rank(datax1,ties.method="first"))
		try(ray0 <- rank(datay0,ties.method="first"))
		try(ray1 <- rank(datay1,ties.method="first"))
			# objects for resorting within split groups
		try(posx <- rank(proby < 0.05,ties.method="first"))
		try(posy <- rank(probx < 0.05,ties.method="first"))
			# objects for resorting the split groups together again
			# F goes first, so combine 0 (p large) before 1 (p small)
		if(length(datax0)>20 & length(datax1)>20){
			try(Locoutx0 <- lapply(as.list(1:length(NTCxN)),function(i){Locfunc2(sort(datax0),modout$center[1],NTCxN[[i]])}))
				# use only neg droplets other channels as predictors
			sigdifx <- mean((quantile((datax1-modout$center[1])[datax1-modout$center[1]<=0])/quantile((datax0-modout$center[1])[datax0-modout$center[1]<=0]))[2:4])
			if(is.na(sigdifx)|(sigdifx<1)){sigdifx <- 1}
			datax1c <- (datax1-modout$center[1])/sigdifx+modout$center[1]
				# difference of variability of the null cluster
			try(Locoutx1 <- lapply(as.list(1:length(NTCxN)),function(i){Locfunc2(sort(datax0),modout$center[1],NTCxN[[i]],sort(datax1c))}))
				# predict for pos droplets other channels after correction

			try(problistx0 <- lapply(Locoutx0,function(x){x$fdr}))
			try(problistx1 <- lapply(Locoutx1,function(x){x$fdr}))
				# store posterior probabilities
			try(datafitx <- mapply(c,lapply(Locoutx0,function(x){x$datafit[rax0]}),lapply(Locoutx1,function(x){x$datafit[rax1]}),SIMPLIFY=F))
			try(datafitx <- lapply(datafitx,function(x){x[posx]}))
				# store data fit sorted as in original data
			try(ntcfitx <- mapply(c,lapply(Locoutx0,function(x){x$datanull[rax0]}),lapply(Locoutx1,function(x){x$datanull[rax1]}),SIMPLIFY=F))
			try(ntcfitx <- lapply(ntcfitx,function(x){x[posx]}))
				# store NTC fit sorted as in original data
			try(pifitsx <- lapply(Locoutx0,function(x){x$pifits}))
				# store pi fits
			mux <- lapply(Locoutx0,function(x){x$res[1]})
			sigx <- lapply(Locoutx0,function(x){x$res[2]})
			pix <- lapply(Locoutx0,function(x){x$res[3]})
Ttime4[[1]] <- Sys.time()
			try(probrx0 <- matrix(unlist(problistx0),ncol=length(NTCxN)))
			try(probmx0 <- rowMeans(probrx0))
			try(probx0 <- BoundedAntiMean(probmx0,w=rep(1,length(datax0)),a=rep(0,length(datax0)),b=rep(1,length(datax0)))[rax0])
				# bounded antitonic regression on the probabilities
				# sorted results are reordered using rax (original ranks)
			try(probrx1 <- matrix(unlist(problistx1),ncol=length(NTCxN)))
			try(probmx1 <- rowMeans(probrx1))
			try(probx1 <- BoundedAntiMean(probmx1,w=rep(1,length(datax1)),a=rep(0,length(datax1)),b=rep(1,length(datax1)))[rax1])
				# bounded antitonic regression on the probabilities
				# sorted results are reordered using rax (original ranks)
			try(probx <- c(probx0,probx1)[posx])
				# partition probabilities
			try(if(length(NTCxN)>1){
				probrx <- rbind(probrx0[rax0,],probrx1[rax1,])[posx,]
				probsx0 <- apply(probrx0,1,sd)
				probciLx0 <- probmx0-qt(0.975,length(NTCxN)-1)*probsx0*sqrt(1+1/length(NTCxN))
				probciUx0 <- probmx0+qt(0.975,length(NTCxN)-1)*probsx0*sqrt(1+1/length(NTCxN))
				probciLx0 <- BoundedAntiMean(probciLx0,w=rep(1,length(datax0)),a=rep(0,length(datax0)),b=rep(1,length(datax0)))[rax0]
				probciUx0 <- BoundedAntiMean(probciUx0,w=rep(1,length(datax0)),a=rep(0,length(datax0)),b=rep(1,length(datax0)))[rax0]
				probsx1 <- apply(probrx1,1,sd)
				probciLx1 <- probmx1-qt(0.975,length(NTCxN)-1)*probsx1*sqrt(1+1/length(NTCxN))
				probciUx1 <- probmx1+qt(0.975,length(NTCxN)-1)*probsx1*sqrt(1+1/length(NTCxN))
				probciLx1 <- BoundedAntiMean(probciLx1,w=rep(1,length(datax1)),a=rep(0,length(datax1)),b=rep(1,length(datax1)))[rax1]
				probciUx1 <- BoundedAntiMean(probciUx1,w=rep(1,length(datax1)),a=rep(0,length(datax1)),b=rep(1,length(datax1)))[rax1]
				probciLx <- c(probciLx0,probciLx1)[posx]
				probciUx <- c(probciUx0,probciUx1)[posx]
			}) # confidence limits, if possible
			try(if(length(NTCxN)==1){
				probrx <- c(probrx0[rax0,],probrx1[rax1,])[posx]})
Ttime4[[2]] <- Sys.time()
		} # end x-loop enough pos/neg data
		else{
			try(datafitx <- lapply(Locoutx,function(x){x$datafit[rax]}))
				# store data fit sorted as in original data
			try(ntcfitx <- lapply(Locoutx,function(x){x$datanull[rax]}))
				# store NTC fit sorted as in original data
			try(pifitsx <- lapply(Locoutx,function(x){x$pifits}))
				# store pi fits
			mux <- lapply(Locoutx,function(x){x$res[1]})
			sigx <- lapply(Locoutx,function(x){x$res[2]})
			pix <- lapply(Locoutx,function(x){x$res[3]})
			try(if(length(NTCxN)>1){
				probrx <- probrx[rax,]
				probsx <- apply(probrx,1,sd)
				probciLx <- probmx-qt(0.975,length(NTCxN)-1)*probsx*sqrt(1+1/length(NTCxN))
				probciUx <- probmx+qt(0.975,length(NTCxN)-1)*probsx*sqrt(1+1/length(NTCxN))
				probciLx <- BoundedAntiMean(probciLx,w=rep(1,length(datax)),a=rep(0,length(datax)),b=rep(1,length(datax)))[rax]
				probciUx <- BoundedAntiMean(probciUx,w=rep(1,length(datax)),a=rep(0,length(datax)),b=rep(1,length(datax)))[rax]
			}) # confidence limits, if possible
			try(if(length(NTCxN)==1){
				probrx <- probrx[rax]})

		} # end x-loop only pos or neg data

		if(length(datay0)>20 & length(datay1)>20){
			try(Locouty0 <- lapply(as.list(1:length(NTCyN)),function(i){Locfunc2(sort(datay0),modout$center[2],NTCyN[[i]])}))
				# use only neg droplets other channels as predictors
			sigdify <- mean((quantile((datay1-modout$center[2])[datay1-modout$center[2]<=0])/quantile((datay0-modout$center[2])[datay0-modout$center[2]<=0]))[2:4])
			if(is.na(sigdify)|(sigdify<1)){sigdify <- 1}
			datay1c <- (datay1-modout$center[2])/sigdify+modout$center[2]
				# difference of variability of the null cluster
			try(Locouty1 <- lapply(as.list(1:length(NTCyN)),function(i){Locfunc2(sort(datay0),modout$center[2],NTCyN[[i]],sort(datay1c))}))
				# predict for pos droplets other channels after correction

			try(problisty0 <- lapply(Locouty0,function(x){x$fdr}))
			try(problisty1 <- lapply(Locouty1,function(x){x$fdr}))
				# store posterior probabilities
			try(datafity <- mapply(c,lapply(Locouty0,function(x){x$datafit[ray0]}),lapply(Locouty1,function(x){x$datafit[ray1]}),SIMPLIFY=F))
			try(datafity <- lapply(datafity,function(x){x[posy]}))
				# store data fit sorted as in original data
			try(ntcfity <- mapply(c,lapply(Locouty0,function(x){x$datanull[ray0]}),lapply(Locouty1,function(x){x$datanull[ray1]}),SIMPLIFY=F))
			try(ntcfity <- lapply(ntcfity,function(x){x[posy]}))
				# store NTC fit sorted as in original data
			try(pifitsy <- lapply(Locouty0,function(x){x$pifits}))
				# store pi fits
			muy <- lapply(Locouty0,function(x){x$res[1]})
			sigy <- lapply(Locouty0,function(x){x$res[2]})
			piy <- lapply(Locouty0,function(x){x$res[3]})
Ttime4[[3]] <- Sys.time()
			try(probry0 <- matrix(unlist(problisty0),ncol=length(NTCyN)))
			try(probmy0 <- rowMeans(probry0))
			try(proby0 <- BoundedAntiMean(probmy0,w=rep(1,length(datay0)),a=rep(0,length(datay0)),b=rep(1,length(datay0)))[ray0])
				# bounded antitonic regression on the probabilities
				# sorted results are reordered using ray (original ranks)
			try(probry1 <- matrix(unlist(problisty1),ncol=length(NTCyN)))
			try(probmy1 <- rowMeans(probry1))
			try(proby1 <- BoundedAntiMean(probmy1,w=rep(1,length(datay1)),a=rep(0,length(datay1)),b=rep(1,length(datay1)))[ray1])
				# bounded antitonic regression on the probabilities
				# sorted results are reordered using ray (original ranks)
			try(proby <- c(proby0,proby1)[posy])
			try(if(length(NTCyN)>1){
				probry <- rbind(probry0[ray0,],probry1[ray1,])[posy,]
				probsy0 <- apply(probry0,1,sd)
				probciLy0 <- probmy0-qt(0.975,length(NTCyN)-1)*probsy0*sqrt(1+1/length(NTCyN))
				probciUy0 <- probmy0+qt(0.975,length(NTCyN)-1)*probsy0*sqrt(1+1/length(NTCyN))
				probciLy0 <- BoundedAntiMean(probciLy0,w=rep(1,length(datay0)),a=rep(0,length(datay0)),b=rep(1,length(datay0)))[ray0]
				probciUy0 <- BoundedAntiMean(probciUy0,w=rep(1,length(datay0)),a=rep(0,length(datay0)),b=rep(1,length(datay0)))[ray0]
				probsy1 <- apply(probry1,1,sd)
				probciLy1 <- probmy1-qt(0.975,length(NTCyN)-1)*probsy1*sqrt(1+1/length(NTCyN))
				probciUy1 <- probmy1+qt(0.975,length(NTCyN)-1)*probsy1*sqrt(1+1/length(NTCyN))
				probciLy1 <- BoundedAntiMean(probciLy1,w=rep(1,length(datay1)),a=rep(0,length(datay1)),b=rep(1,length(datay1)))[ray1]
				probciUy1 <- BoundedAntiMean(probciUy1,w=rep(1,length(datay1)),a=rep(0,length(datay1)),b=rep(1,length(datay1)))[ray1]
				probciLy <- c(probciLy0,probciLy1)[posy]
				probciUy <- c(probciUy0,probciUy1)[posy]
			}) # confidence limits, if possible
			try(if(length(NTCyN)==1){
				probry <- c(probry0[ray0,],probry1[ray1,])[posy]
})

Ttime4[[4]] <- Sys.time()
		} # end y-loop enough pos/neg data
		else{
			try(datafity <- lapply(Locouty,function(x){x$datafit[ray]}))
				# store data fit sorted as in original data
			try(ntcfity <- lapply(Locouty,function(x){x$datanull[ray]}))
				# store NTC fit sorted as in original data
			try(pifitsy <- lapply(Locouty,function(x){x$pifits}))
				# store pi fits
			muy <- lapply(Locouty,function(x){x$res[1]})
			sigy <- lapply(Locouty,function(x){x$res[2]})
			piy <- lapply(Locouty,function(x){x$res[3]})
			try(if(length(NTCyN)>1){
				probry <- probry[ray,]
				probsy <- apply(probry,1,sd)
				probciLy <- probmy-qt(0.975,length(NTCyN)-1)*probsy*sqrt(1+1/length(NTCyN))
				probciUy <- probmy+qt(0.975,length(NTCyN)-1)*probsy*sqrt(1+1/length(NTCyN))
				probciLy <- BoundedAntiMean(probciLy,w=rep(1,length(datay)),a=rep(0,length(datay)),b=rep(1,length(datay)))[ray]
				probciUy <- BoundedAntiMean(probciUy,w=rep(1,length(datay)),a=rep(0,length(datay)),b=rep(1,length(datay)))[ray]
			}) # confidence limits, if possible
			try(if(length(NTCyN)==1){
				probry <- probry[ray]})
		} # end y-loop only pos or neg data
	} # end NTC-version loop

# non-NTC version fit
	else{
		# the same as NTC version but with Locfunc1
		try(datax0 <- datax[proby >= 0.8])
		try(datax1 <- datax[proby < 0.8])
		try(datay0 <- datay[probx >= 0.8])
		try(datay1 <- datay[probx < 0.8])
			# data split over other channel
		try(rax0 <- rank(datax0,ties.method="first"))
		try(rax1 <- rank(datax1,ties.method="first"))
		try(ray0 <- rank(datay0,ties.method="first"))
		try(ray1 <- rank(datay1,ties.method="first"))
			# objects for resorting within split groups
		try(posx <- rank(proby < 0.8,ties.method="first"))
		try(posy <- rank(probx < 0.8,ties.method="first"))
			# locx
		if(length(datax0)>20 & length(datax1)>20){
			try(Locoutx0 <- Locfunc1(sort(datax0),mux))
				# use only neg droplets other channels as predictors
			sigdifx <- mean((quantile((datax1-modout$center[1])[datax1-modout$center[1]<=0])/quantile((datax0-modout$center[1])[datax0-modout$center[1]<=0]))[2:4])
			if(is.na(sigdifx)|(sigdifx<1)){sigdifx <- 1}
			datax1c <- (datax1-modout$center[1])/sigdifx+modout$center[1]
				# difference of variability of the null cluster
			try(Locoutx1 <- Locfunc1(sort(datax0),modout$center[1],sort(datax1c)))
				# predict for pos droplets other channels after correction
			try(mux <- Locoutx0$res[1])
			try(sigx <- Locoutx0$res[2])
			try(pix <- Locoutx0$res[3])
			try(probrx0 <- Locoutx0$fdr)
			try(probrx1 <- Locoutx1$fdr)
				# store posterior probabilities
			try(datafitx <- list(datafit=c(Locoutx0$datafit[rax0],Locoutx1$datafit[rax1])[posx]))
				# store datafit sorted as in original data
			try(ntcfitx <- list(nullfit=c(Locoutx0$datanull[rax0],Locoutx1$datanull[rax1])[posx]))
				# store NTC fit sorted as in original data
			try(pifitsx <- Locoutx0$pifits)
Ttime4[[1]] <- Sys.time()
			try(probx0 <- BoundedAntiMean(probrx0,w=rep(1,length(datax0)),a=rep(0,length(datax0)),b=rep(1,length(datax0)))[rax0])
			try(probx1 <- BoundedAntiMean(probrx1,w=rep(1,length(datax1)),a=rep(0,length(datax1)),b=rep(1,length(datax1)))[rax1])
				# bounded antitonic regression on the probabilities
				# sorted results are reordered using rax (original ranks)
			try(probx <- c(probx0,probx1)[posx])
			try(probrx <- c(probrx0[rax0],probrx1[rax1])[posx])
				# partition probabilities
Ttime4[[2]] <- Sys.time()
		} else{
			try(datafitx <- list(datafit=Locoutx$datafit[rax]))
				# store data fit sorted as in original data
			try(ntcfitx <- list(nullfit=Locoutx$datanull[rax]))
				# store NTC fit sorted as in original data
			try(pifitsx <- Locoutx$pifits)
				# store pi fits
			mux <- Locoutx$res[1]
			sigx <- Locoutx$res[2]
			pix <- Locoutx$res[3]
		} # end y-loop only pos or neg data
			# locy
		if(length(datay0)>20 & length(datay1)>20){
			try(Locouty0 <- Locfunc1(sort(datay0),muy))
				# use only neg droplets other channels as predictors
			sigdify <- mean((quantile((datay1-modout$center[2])[datay1-modout$center[2]<=0])/quantile((datay0-modout$center[2])[datay0-modout$center[2]<=0]))[2:4])
			if(is.na(sigdify)|(sigdify<1)){sigdify <- 1}
			datay1c <- (datay1-modout$center[2])/sigdify+modout$center[2]
				# difference of variability of the null cluster
			try(Locouty1 <- Locfunc1(sort(datay0),modout$center[2],sort(datay1c)))
				# predict for pos droplets other channels after correction
			try(muy <- Locouty0$res[1])
			try(sigy <- Locouty0$res[2])
			try(piy <- Locouty0$res[3])
			try(probry0 <- Locouty0$fdr)
			try(probry1 <- Locouty1$fdr)
				# store posterior probabilities
			try(datafity <- list(datafit=c(Locouty0$datafit[ray0],Locouty1$datafit[ray1])[posy]))
				# store datafit sorted as in original data
			try(ntcfity <- list(nullfit=c(Locouty0$datanull[ray0],Locouty1$datanull[ray1])[posy]))
				# store NTC fit sorted as in original data
			try(pifitsy <- Locouty0$pifits)
Ttime4[[1]] <- Sys.time()
			try(proby0 <- BoundedAntiMean(probry0,w=rep(1,length(datay0)),a=rep(0,length(datay0)),b=rep(1,length(datay0)))[ray0])
			try(proby1 <- BoundedAntiMean(probry1,w=rep(1,length(datay1)),a=rep(0,length(datay1)),b=rep(1,length(datay1)))[ray1])
				# bounded antitonic regression on the probabilities
				# sorted results are reordered using ray (original ranks)
			try(proby <- c(proby0,proby1)[posy])
			try(probry <- c(probry0[ray0],probry1[ray1])[posy])
				# partition probabilities
Ttime4[[4]] <- Sys.time()
	}	else{
			try(datafity <- list(datafit=Locouty$datafit[ray]))
				# store data fit sorted as in original data
			try(ntcfity <- list(nullfit=Locouty$datanull[ray]))
				# store NTC fit sorted as in original data
			try(pifitsy <- Locouty$pifits)
				# store pi fits
			muy <- Locouty$res[1]
			sigy <- Locouty$res[2]
			piy <- Locouty$res[3]
		} # end y-loop only pos or neg data
	}

# write output
	try(prl <- list(mu=cbind(unlist(mux),unlist(muy)),sig=cbind(unlist(sigx),unlist(sigy)),pi=cbind(unlist(pix),unlist(piy)),sigdif=c(sigdifx,sigdify)))
	try(px <- min(mean(unlist(pix)),1))
	try(py <- min(mean(unlist(piy)),1))
	try(prob2d <- cbind(probx,proby))
		# posterior probabilities to be negative per well
	piz <- NULL
	try(piz <- cbind(probx*proby,probx*(1-proby),(1-probx)*proby,(1-probx)*(1-proby)))
		# posterior probabilities per cluster
	ndropx <- length(probx)
	nposx <- sum(probx<0.05)
	nnegx <- sum(probx>0.8)
	nrainx <- ndropx-nposx-nnegx
	ndropy <- length(proby)
	nposy <- sum(proby<0.05)
	nnegy <- sum(proby>0.8)
	nrainy <- ndropy-nposy-nnegy
	pposx <- nposx/ndropx
	pposy <- nposy/ndropy

	if(length(NTCxN)>1){
		psdx <- sd(pmin(unlist(pix),1))
		csdx <- psdx/px  # delta rule with g'(p) = 1/p^2
		if(px==0 & psdx==0){csdx <- 0}
		pircix <- pmin(pmax(mean(unlist(pix))+c(-1,1)*qnorm(0.975)*psdx*sqrt(1+1/length(NTCxN)),0),1)
		concix <- pmax(-log(px)+c(-1,1)*qnorm(0.975)*sqrt(csdx^2*(1+1/length(NTCxN))+(1/px-1)/ndropx),0)
		posrx <- apply(probrx,2,function(x){sum(x<0.05)})
		negrx <- apply(probrx,2,function(x){sum(x>0.8)})
		rainrx <- ndropx-posrx-negrx
		ciposx <- round(pmin(pmax(mean(posrx)+c(-1,1)*qnorm(0.975)*sd(posrx),0),ndropx))
		cinegx <- round(pmin(pmax(mean(negrx)+c(-1,1)*qnorm(0.975)*sd(negrx),0),ndropx))
		cirainx <- round(pmin(pmax(mean(rainrx)+c(-1,1)*qnorm(0.975)*sd(rainrx),0),ndropx))
		if(cirainx[1] > nrainx){cirainx[1] <- nrainx}
		if(cirainx[2] < nrainx){cirainx[2] <- nrainx}
		psdposx <- sd(posrx/ndropx)
		csdposx <- psdposx/(1-pposx)
		if(pposx==0 & psdposx==0){csdposx <- 0}
		concposx <- pmax(-log(1-pposx)+c(-1,1)*qnorm(0.975)*sqrt(csdposx^2*(1+1/length(NTCxN))+(1/(1-pposx)-1)/ndropx),0)
	}
	else{ pircix <- c(0,1)
		concix <- c(0,Inf)
		ciposx <- c(0,ndropx)
		cinegx <- c(0,ndropx)
		cirainx <- c(0,ndropx)
		psdx <- NA
		csdx <- NA
		psdposx <- NA
		csdposx <- NA
		concposx <- c(0,Inf)}

	if(length(NTCyN)>1){
		psdy <- sd(pmin(unlist(piy),1))
		csdy <- psdy/py  # delta rule with g'(p) = 1/p^2
		if(py==0 & psdy==0){csdy <- 0}
		pirciy <- pmin(pmax(mean(unlist(piy))+c(-1,1)*qnorm(0.975)*psdy*sqrt(1+1/length(NTCyN)),0),1)
		conciy <- pmax(-log(py)+c(-1,1)*qnorm(0.975)*sqrt(csdy^2*(1+1/length(NTCyN))+(1/py-1)/ndropy),0)
		posry <- apply(probry,2,function(y){sum(y<0.05)})
		negry <- apply(probry,2,function(y){sum(y>0.8)})
		rainry <- ndropy-posry-negry
		ciposy <- round(pmin(pmax(mean(posry)+c(-1,1)*qnorm(0.975)*sd(posry),0),ndropy))
		cinegy <- round(pmin(pmax(mean(negry)+c(-1,1)*qnorm(0.975)*sd(negry),0),ndropy))
		cirainy <- round(pmin(pmax(mean(rainry)+c(-1,1)*qnorm(0.975)*sd(rainry),0),ndropy))
		if(cirainy[1] > nrainy){cirainy[1] <- nrainy}
		if(cirainy[2] < nrainy){cirainy[2] <- nrainy}
		pposy <- nposy/ndropy
		psdposy <- sd(posry/ndropy)
		csdposy <- psdposy/(1-pposy)
		if(pposy==0 & psdposy==0){csdposy <- 0}
		concposy <- pmax(-log(1-pposy)+c(-1,1)*qnorm(0.975)*sqrt(csdposy^2*(1+1/length(NTCyN))+(1/(1-pposy)-1)/ndropy),0)
	}
	else{ pirciy <- c(0,1)
		conciy <- c(0,Inf)
		ciposy <- c(0,ndropy)
		cinegy <- c(0,ndropy)
		cirainy <- c(0,ndropy)
		psdy <- NA
		csdy <- NA
		psdposy <- NA
		csdposy <- NA
		concposy <- c(0,Inf)}

	try(concx <- rbind(c(1-px,-log(px)*1000/vol,pposx,-log(1-pposx)*1000/vol),
				c(psdx,csdx*1000/vol,psdposx,csdposx*1000/vol),
				c(1-pircix[1],concix[1]*1000/vol,ciposx[1]/ndropx,concposx[1]*1000/vol),
				c(1-pircix[2],concix[2]*1000/vol,ciposx[2]/ndropx,concposx[2]*1000/vol)))
	try(colnames(concx) <- c("prob_robust","conc_robust","prob_thres","conc_thres"))
	try(rownames(concx) <- c("est","sd","CI LB","CI UB"))
		# column 1: probability to be positive with robust estimator
		# column 2: concentration with robust estimator
		# column 3: probability to be positive with 1% threshold
		# column 4: concentration with 1% threshold
		# row 1: estimators
		# row 2: standard deviations
		# row 3: lower bound confidence intervals (including Poisson)
		# row 4: upper bound confidence intervals (including Poisson)

	try(treshx <- rbind(c(nposx/ndropx,ciposx/ndropx,nposx,ciposx),
				c(nnegx/ndropx,cinegx/ndropx,nnegx,cinegx),
				c(nrainx/ndropx,cirainx/ndropx,nrainx,cirainx)))
	try(colnames(treshx) <- c("prob","CI low p","CI high p","count","CI low n","CI high n"))
	try(rownames(treshx) <- c("pos","neg","rain"))

	try(concy <- rbind(c(1-py,-log(py)*1000/vol,pposy,-log(1-pposy)*1000/vol),
				c(psdy,csdy*1000/vol,psdposy,csdposy*1000/vol),
				c(1-pirciy[1],conciy[1]*1000/vol,ciposy[1]/ndropy,concposy[1]*1000/vol),
				c(1-pirciy[2],conciy[2]*1000/vol,ciposy[2]/ndropy,concposy[2]*1000/vol)))
	try(colnames(concy) <- c("prob_robust","conc_robust","prob_thres","conc_thres"))
	try(rownames(concy) <- c("est","sd","CI LB","CI UB"))
		# overall probability to be / proportion negative
		# concentration
	try(treshy <- rbind(c(nposy/ndropy,ciposy/ndropy,nposy,ciposy),
				c(nnegy/ndropy,cinegy/ndropy,nnegy,cinegy),
				c(nrainy/ndropy,cirainy/ndropy,nrainy,cirainy)))
	try(colnames(treshy) <- c("prob","CI low p","CI high p","count","CI low n","CI high n"))
	try(rownames(treshy) <- c("pos","neg","rain"))



# convergence message
	time2 <- Sys.time()
	timerun <- time2-time1

# cat("Timetest\n")
# cat("Initial objects:",Ttime1-time1,"\n")
# cat("Initial locfunc:",Ttime2-Ttime1,"\n")
# cat("Initial isotonic:",Ttime3-Ttime2,"\n")
# cat("x channel locfunc PN:",Ttime4[[1]]-Ttime3,"\n")
# cat("x channel isotonic PN:",Ttime4[[2]]-Ttime4[[1]],"\n")
# cat("y channel locfunc PN:",Ttime4[[3]]-Ttime4[[2]],"\n")
# cat("y channel isotonic PN:",Ttime4[[4]]-Ttime4[[3]],"\n")
# cat("write output:",time2-Ttime4[[4]],"\n")

	cat("Estimation",modout$name,"completed in",timerun,attr(timerun,"units"),". \n")

# Print the results

cat("\n")
cat("Well",modout$name,":\n")
cat("\n")
cat("  Channel 1:\n")
cat(ndropx,"partitions\n")
cat("Robust estimator:\n")
cat("  ",round(ndropx*concx[1,1]),"positive, p =",round(concx[1,1],4),"( CI: [",round(concx[3,1],4),",",round(concx[4,1],4),"] )\n")
cat("  ",round(concx[1,2],1),"copies/mul ( CI: [",round(concx[3,2],1),",",round(concx[4,2],1),"] )\n")
cat("Threshold estimator:\n")
cat("  ",round(ndropx*concx[1,3]),"positive, p =",round(concx[1,3],4),"( CI: [",round(concx[3,3],4),",",round(concx[4,3],4),"] )\n")
cat("  ",round(concx[1,4],1),"copies/mul ( CI: [",round(concx[3,4],1),",",round(concx[4,4],1),"] )\n")
cat("\n")
cat("  Channel 2:\n")
cat(ndropy,"partitions\n")
cat("Robust estimator:\n")
cat("  ",round(ndropy*concy[1,1]),"positive, p =",round(concy[1,1],4),"( CI: [",round(concy[3,1],4),",",round(concy[4,1],4),"] )\n")
cat("  ",round(concy[1,2],1),"copies/mul ( CI: [",round(concy[3,2],1),",",round(concy[4,2],1),"] )\n")
cat("Threshold estimator:\n")
cat("  ",round(ndropy*concy[1,3]),"positive, p =",round(concy[1,3],4),"( CI: [",round(concy[3,3],4),",",round(concy[4,3],4),"] )\n")
cat("  ",round(concy[1,4],1),"copies/mul ( CI: [",round(concy[3,4],1),",",round(concy[4,4],1),"] )\n")
cat("\n")
cat("CI's for concentration include Poisson variability.\n")
cat("\n")
cat("__________________________________________________\n")

rlist <- list(conc=list(x=concx,y=concy),tresh=list(x=treshx,y=treshy),fits=prl,data=data2d,droppi=prob2d,dropci=list(x=cbind(probciLx,probciUx),y=cbind(probciLy,probciUy)),reppi=list(x=probrx,y=probry),densfits=list(x=datafitx,y=datafity),ntcfits=list(x=ntcfitx,y=ntcfity),pifits=list(x=pifitsx,y=pifitsy),cluspi=piz)
return(rlist)
}



#######################
######   Shell   ######
#######################

if(any(duplicated(names(datalist)))){
stop("Duplicated names in data, please provide unique names for each reaction")}

### Preprocessing

# Identify 3 important clusters

nwell <- length(datalist)
dataALL <- NULL
for(i in 1:nwell){
newdata <- data.frame(cbind(datalist[[i]],i))
names(newdata) <- c(names(datalist[[1]]),"id")
dataALL <- rbind(dataALL,newdata)}
	# make a dataframe with all the data (being careful with possible names)

nullclus <- NULL
	# nullclus will be a vector with the id's of 3 clusters (00,10,01)
clinm <- NULL
	# clinm is the method used to determine initial clusters

	# option for prior clusters given in the data
if(max(names(dataALL)=="Cluster")==1){
	if(nlevels(factor(dataALL$Cluster))==4){
		x1t <- dataALL[dataALL$Cluster==1,]
		x2t <- dataALL[dataALL$Cluster==2,]
		x3t <- dataALL[dataALL$Cluster==3,]
		x4t <- dataALL[dataALL$Cluster==4,]
		clus <- floor(rank(c(median(x1t[,1]),median(x2t[,1]),median(x3t[,1]),median(x4t[,1])))/3)+
		2*floor(rank(c(median(x1t[,2]),median(x2t[,2]),median(x3t[,2]),median(x4t[,2])))/3)+1
		nullclus <- c(which(clus==1),which(clus==2),which(clus==3))
		clinm <- "input"
		
		# check NTC data. If NTC are not in first cluster, something is wrong with prior clusters.
		clusNTC <- NULL
		if(!is.null(NTC)){
			for(i in 1:length(NTC)){
				if(any(names(datalist)==NTC[i])){
					NTCclus <- datalist[[which(names(datalist)==NTC[i])]]$Cluster
					if(mean(NTCclus==nullclus[1])<0.9){clinm <- NULL
						cat("      Warning: NTC data not in prior double negative cluster for",NTC[i],"\n")}
	}}}}
	if(is.null(clinm)){nullclus <- NULL
		cat("      Note: prior clusters will not be used.\n")}
	if(length(nullclus)!=3){nullclus <- NULL
		cat("      Note: prior clusters will not be used.\n")}
}

	# when initial clusters are not used, define somewhat arbitrary cutpoints
if(is.null(nullclus)){
	cutx=mean(quantile(dataALL[,1],c(0.01,0.999)))
	cuty=mean(quantile(dataALL[,2],c(0.01,0.999)))
		# arbitrary cutpoint
	clinm <- "cutori"
	if(!is.null(NTC)){
		for(i in 1:length(NTC)){
			if(any(names(datalist)==NTC[i])){
					NTCdat <- datalist[[which(names(datalist)==NTC[i])]]
					cutn <- median(NTCdat[,1])+5*huber(NTCdat[,1])$s
					if(cutx<cutn){cutx <- cutn
						clinm <- "cutcor"}
					cutn <- median(NTCdat[,2])+5*huber(NTCdat[,2])$s
					if(cuty<cutn){cuty <- cutn
						clinm <- "cutcor"}
	}}}
		# make sure NTC's are all below cutpoint
	Clust <- 1+(dataALL[,1]>cutx)+(dataALL[,2]>cuty)*2
	nullclus <- c(1,2,3)
	for(i in 1:length(names(datalist))){
		datalist[[i]] <- data.frame(datalist[[i]][,1:2],Cluster=Clust[dataALL$id==i])}
}
names(nullclus) <- c("clus00","clus10","clus01")


# Identify cluster centres

modout <- lapply(datalist,modeRetFin)
names(modout) <- names(datalist)
for(i in 1:length(datalist)){modout[[i]]$name <- names(datalist)[[i]]}

# separate procedure for NTC
# centre double negative must be mode, other centres must be NA
if(!is.null(NTC)){
	for(i in 1:length(NTC)){
		if(any(names(modout)==NTC[i])){
			modout[[which(names(modout)==NTC[i])]]$diffc <- F
			modout[[which(names(modout)==NTC[i])]]$center[3:6] <- NA
			NTCxcent <- mlv(datalist[[which(names(datalist)==NTC[i])]][,1],method="asselin")$M
			NTCycent <- mlv(datalist[[which(names(datalist)==NTC[i])]][,2],method="asselin")$M
			modout[[which(names(modout)==NTC[i])]]$center[1:2] <- c(NTCxcent,NTCycent)
}}}

diffC <- sapply(modout,function(x){x$diffc})
if(sum(diffC)>0){
cat("      Possible differences found between well methods and plate methods in the following wells:\n")
cat("      ",names(modout)[diffC],"\n")
}
if(mean(diffC)>0.1 & Plate==T & clinm!="input"){
cat("      Use option Plate=F if fluorescence output is not consistent over all wells in dataset.\n")}



# Checking and correcting cluster centres

	# all univariate double negative centres in a vector
C1x <- sapply(modout,function(x){x$center[1]})
C1y <- sapply(modout,function(x){x$center[2]})
	# all univariate distances single positive to double negative centres in a vector
C2x <- sapply(modout,function(x){x$center[3]-x$center[1]})
C2y <- sapply(modout,function(x){x$center[4]-x$center[2]})
C3x <- sapply(modout,function(x){x$center[5]-x$center[1]})
C3y <- sapply(modout,function(x){x$center[6]-x$center[2]})
	# medians and robust standard deviations of all non-NA centres/distances
Cmat <- cbind(rep(NA,6),rep(NA,6))
if(sum(!is.na(C1x))>1){Cmat[1,] <- c(median(C1x,na.rm=T),huber(C1x[!is.na(C1x)])$s)}
if(sum(!is.na(C1y))>1){Cmat[2,] <- c(median(C1y,na.rm=T),huber(C1y[!is.na(C1y)])$s)}
if(sum(!is.na(C2x))>1){Cmat[3,] <- c(median(C2x,na.rm=T),huber(C2x[!is.na(C2x)])$s)}
if(sum(!is.na(C2y))>1){Cmat[4,] <- c(median(C2y,na.rm=T),huber(C2y[!is.na(C2y)])$s)}
if(sum(!is.na(C3x))>1){Cmat[5,] <- c(median(C3x,na.rm=T),huber(C3x[!is.na(C3x)])$s)}
if(sum(!is.na(C3y))>1){Cmat[6,] <- c(median(C3y,na.rm=T),huber(C3y[!is.na(C3y)])$s)}

	# correct only if user allows for using plate info
if((Plate==T) & (length(datalist) > 2)){modplat <- lapply(modout,checkmod1)}
else{modplat <- modout}
names(modplat) <- names(modout)

	# impute NTC single positive centres for rotation step.
if(!is.null(NTC)){
	for(i in 1:length(NTC)){
		if(any(names(modplat)==NTC[i])){
			modplat[[which(names(modplat)==NTC[i])]]$center[3:6] <- rep(modplat[[which(names(modplat)==NTC[i])]]$center[1:2],2)+Cmat[3:6,1]
	}}}

if((max(is.na(Cmat))==1) & ((Rotate[1]=="Plate")|(Rotate[1]==T))){Rotate<-F
cat("      Warning: rotation not possible.")}

	# execute the actual rotation
if(Rotate[1]=="Plate"){modplatA <- modplat
	for(i in 1:nwell){modplat[[i]]$center[3:6] <- rep(modplatA[[i]]$center[1:2],2)+Cmat[3:6,1]}
	datarotlist <- lapply(as.list(1:length(datalist)),function(i){Rotfunc(datalist[[i]],modplat[[i]])})
} else if(any(names(modout)%in%Rotate))
	{modplatB <- modplat
	refs <- which(names(modout)%in%Rotate)
	Cw <- NULL
	for(i in 1:length(refs)){
	Cw[[i]] <- modout[[refs[i]]]$center[3:6]-rep(modplat[[refs[i]]]$center[1:2],2)}
	Cwells <- matrix(rowMeans(matrix(unlist(Cw),ncol=length(refs))),nrow=2)
	for(i in 1:nwell){modplat[[i]]$center[3:6] <- rep(modplatB[[i]]$center[1:2],2)+Cwells}
	datarotlist <- lapply(as.list(1:length(datalist)),function(i){Rotfunc(datalist[[i]],modplat[[i]])})
} else if(Rotate==T){datarotlist <- lapply(as.list(1:length(datalist)),function(i){Rotfunc(datalist[[i]],modplat[[i]])})
} else {datarotlist <- datalist
	cat("      Warning: not checking for rotation can lead to wrong results if rotation is present.\n")}

names(datarotlist) <- names(datalist)

cat("Preprocessing complete\n")

	# create univariate NTC objects, all NTC's in one object
NTCx <- NULL
NTCy <- NULL
NTCxN <- NULL
NTCyN <- NULL
if(!is.null(NTC)){
	NTCx <- list()
	NTCy <- list()
	for(i in 1:length(NTC)){
		if(any(names(datarotlist)==NTC[i])){
			NTCx <- append(NTCx,list(datarotlist[[which(names(datarotlist)==NTC[i])]][,1]-modplat[[which(names(modplat)==NTC[i])]]$center[1]))
			NTCy <- append(NTCy,list(datarotlist[[which(names(datarotlist)==NTC[i])]][,2]-modplat[[which(names(modplat)==NTC[i])]]$center[2]))
				# NTC fluorescence minus null mode (double negative centre)
	}}
		# if there are NTC's, the NTCx object includes all NTC channel data from the first channel
		# if there are NTC's, the NTCy object includes all NTC channel data from the second channel
	if((!is.null(NTCx))&(!is.null(NTCy))){
		NTCxN <- lapply(NTCx,NTCfunc)
		NTCyN <- lapply(NTCy,NTCfunc)
			# generate NTC property data
		cat("NTC identification complete\n")}
	else{cat("      Warning: no wells found with given NTC names.\n")
		cat("NTC identification skipped\n")}}
else{cat("NTC identification skipped\n")}

	# This is the real procedure, all calculations happen internally in this step
resultlist <- lapply(as.list(1:length(datarotlist)),function(i){Shellfunc(datarotlist[[i]],modplat[[i]])})
names(resultlist) <- names(datarotlist)
for(i in 1:length(datarotlist)){resultlist[[i]]$name <- names(datarotlist)[[i]]
rawdata <- list(datalist[[i]][,1:2])
names(rawdata) <- "rawdata"
resultlist[[i]] <- append(resultlist[[i]],rawdata)}

cat("Estimation complete\n")

return(resultlist)
}

