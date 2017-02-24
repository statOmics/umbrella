PartClass1d <- function(datalist,vol=0.85,NTC=NULL){
	# datalist must be a list of data: dataframes or vectors.
	# 	Each dataframe has 1 channel of fluorescence intensities
	# vol is the volume of one droplet / partition in nL (0.89 for Bio-Rad's QX-100)
	# NTC is a vector of names indicating which elements from datalist are NTC's
	# 	The names in NTC must match the names of items in the datalist.
require(MASS)
require(mgcv)
require(robust)
require(modeest)
require(OrdMonReg)

########################
######   Proc 1   ######
########################

### function to find the cluster centres

modeRet=function(x,NTC=FALSE,prob=1e-3,method="asselin")
{
if (NTC) return(mlv(x,method=method)$M) 
#mvl implements several methods to find the mode of a univariate unimodal distribution (and thus the most abundant mode of a univariate multimodale distributie)

else 
{
#use cut point to find two modi
cut=mean(quantile(x,c(prob,1-prob)))
hlp=c(mlv(x,method=method)$M,mlv(x[x<cut],method=method)$M,mlv(x[x>cut],method=method)$M,NA)
#minimal mode is null mode
if (hlp[1]/mean(hlp[1:3])>1.1) hlp[4]=min(hlp[1:3]) else hlp[4]=hlp[1]
return(list(mode=hlp[4],name=NULL))
}
}

########################
######   Proc 2   ######
########################

### function to check null modes based on plate information

checkmod1 <- function(modsds){
	modnew <- modsds
		# give warning when a cluster (centre) is not defined.
		# don't warn if it's the single positive of an NTC.
	if(is.na(modsds$mode)){
		modnew$mode <- platmod[1]
		cat("      Warning: null mode in",modsds$name,"not defined.  Plate info used instead.\n")}
	dists <- abs(modsds$mode-platmod[1])/platmod[2]
		# relative difference from plate median (number of robust sd's)
	dists[is.na(dists)] <- 0
	if(dists>5)
		{cat("      Warning: null mode in",modsds$name,"may deviate from other modes in plate. Please check input.\n")}
	return(modnew)
		}

########################
######   Proc 3   ######
########################

# Create NTC object with information

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
######   Proc 4   ######
########################

# Fit spline and normal null distribution (2nd order polynomial) on the data
# back-up procedure for no NTC's

Locfunc1 <- function(datafit,mu){
	# data should be a datavector
	# mu should be the null mode

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
f <- predict(fmod,newdata=data.frame(spmids=datafit),type="response")

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

try(nullfit <- dnorm(datafit,mufit+mu,sigfit)*pifit)
try(nullfit[datafit < mufit+mu-sigfit] <- datasplf[datafit < mufit+mu-sigfit])
res <- c(mufit+mu,sigfit,pifit)
try(fdr <- nullfit/datasplf)
return(list(fdr=fdr,datanull=nullfit,datafit=datasplf,pifits=pifit,res=res))
}

########################
######   Proc 5   ######
########################

# Fit NTC on the data

Locfunc2 <- function(datafit,mu,NTC){
	# data should be a datavector
	# mu should be the null mode
	# NTC should be a list with NTC model:
	# NTC$NTC are the NTC frequencies
	# NTC$NTCden is the correction to get a density
	# NTC$NTCsig is the robust deviation (mad) of the NTC

	fdr <- NULL
	fmarg <- NULL
	f0 <- NULL
	pi0 <- NULL
	prob <- NULL
	sigma <- NTC$NTCsig
### spline fit

try(if(length(datafit)>3){
	data <- datafit-mu
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

	fmarg <- predict(fitHlp,data.frame(mids=data,dummy=1),type="response")
	f0 <- prob*predict(fitHlp,data.frame(mids=data,dummy=0),type="response")
	fdr <- f0/fmarg
	fdr[data < -sigma] <- 1
}
else{cat("      Warning: not enough data to fit null distribution\n")})

res <- c(mu,NTC$NTCsig,prob)
return(list(fdr=fdr,datafit=fmarg,datanull=f0,pifits=pi0,res=res))
}

########################
######   Proc 6   ######
########################

# Shell to run proc 4, 5

Shellfunc <- function(data,modout){
	time1 <- Sys.time()
	Locoutx <- list()
		# create empty objects to store results
	datax <- sort(data)
	rax <- rank(data,ties.method="first")
	orx <- order(data)
		# univariate channel data and sort objects
	probx <- rep(NA,length(datax))
	problist <- list(x=probx)
	probciL <- probx
	probciU <- probx
	probr <- probx
		# create empty objects to store results

Ttime1 <- Sys.time()

# locfit
	if(!is.null(NTCxN)){
		try(Locoutx <- lapply(as.list(1:length(NTCxN)),function(i){Locfunc2(datax,modout$mode,NTCxN[[i]])}))
			# if NTC's are present, use Locfunc2 for fitting
		try(problist <- lapply(Locoutx,function(x){x$fdr}))
			# store posterior probabilities
		try(datafit <- lapply(Locoutx,function(x){x$datafit[rax]}))
			# store data fit sorted as in original data
		try(ntcfit <- lapply(Locoutx,function(x){x$datanull[rax]}))
			# store NTC fit sorted as in original data
		try(pifits <- lapply(Locoutx,function(x){x$pifits}))
			# store pi fits
		mux <- lapply(Locoutx,function(x){x$res[1]})
		sigx <- lapply(Locoutx,function(x){x$res[2]})
		pix <- lapply(Locoutx,function(x){x$res[3]})
	Ttime2 <- Sys.time()
		try(probr <- matrix(unlist(problist),ncol=length(NTCxN)))
		try(probm <- rowMeans(probr))
		try(probx <- BoundedAntiMean(probm,w=rep(1,length(datax)),a=rep(0,length(datax)),b=rep(1,length(datax)))[rax])
			# bounded antitonic regression on the probabilities
			# sorted results are reordered using rax (original ranks)
		try(if(length(NTCxN)>1){
			probs <- apply(probr,1,sd)
			probciL <- probm-qt(0.975,length(NTCxN)-1)*probs*sqrt(1+1/length(NTCxN))
			probciU <- probm+qt(0.975,length(NTCxN)-1)*probs*sqrt(1+1/length(NTCxN))
			probciL <- BoundedAntiMean(probciL,w=rep(1,length(datax)),a=rep(0,length(datax)),b=rep(1,length(datax)))[rax]
			probciU <- BoundedAntiMean(probciU,w=rep(1,length(datax)),a=rep(0,length(datax)),b=rep(1,length(datax)))[rax]
		})
		try(probr <- probr[rax,])
		}
	else{ try(Locoutx <- Locfunc1(datax,modout$mode))
			# if no NTC's are present, use Locfunc1 for fitting
		try(probr <- Locoutx$fdr)
			# store posterior probabilities
		mux <- Locoutx$res[1]
		sigx <- Locoutx$res[2]
		pix <- Locoutx$res[3]
		try(datafit <- Locoutx$datafit)
			# store data fit sorted as in original data
		try(ntcfit <- Locoutx$datanull)
			# store NTC fit sorted as in original data
	Ttime2 <- Sys.time()
		try(probx <- BoundedAntiMean(probr,rep(1,length(datax)),a=rep(0,length(datax)),b=rep(1,length(datax)))[rax])
			# bounded antitonic regression on the probabilities
			# sorted results are reordered using rax (original ranks)
		try(probr <- probr[rax])
		}

Ttime3 <- Sys.time()

# write output
	try(prl <- list(mu=unlist(mux),sig=unlist(sigx),pi=unlist(pix)))
		# internal estimates
	px <- min(mean(prl$pi),1)
		# mean null probability

	ndrop <- length(probx)
	npos <- sum(probx<0.05)
	nneg <- sum(probx>0.8)
	nrain <- ndrop-npos-nneg
	ppos <- npos/ndrop

	if(length(NTCxN)>1){
		psd <- sd(pmin(prl$pi,1))
		csd <- psd/px  # delta rule with g'(p) = 1/p^2
		if(mean(prl$pi)==0 & psd==0){csd <- 0}
		pirci <- pmin(pmax(mean(prl$pi)+c(-1,1)*qnorm(0.975)*psd*sqrt(1+1/length(NTCxN)),0),1)
		conci <- pmax(-log(px)+c(-1,1)*qnorm(0.975)*sqrt(csd^2*(1+1/length(NTCxN))+(1/px-1)/ndrop),0)
		posr <- apply(probr,2,function(x){sum(x<0.05)})
		negr <- apply(probr,2,function(x){sum(x>0.8)})
		rainr <- ndrop-posr-negr
		cipos <- round(pmin(pmax(mean(posr)+c(-1,1)*qnorm(0.975)*sd(posr),0),ndrop))
		cineg <- round(pmin(pmax(mean(negr)+c(-1,1)*qnorm(0.975)*sd(negr),0),ndrop))
		cirain <- round(pmin(pmax(mean(rainr)+c(-1,1)*qnorm(0.975)*sd(rainr),0),ndrop))
		if(cirain[1] > nrain){cirain[1] <- nrain}
		if(cirain[2] < nrain){cirain[2] <- nrain}
		psdpos <- sd(posr/ndrop)
		csdpos <- psdpos/(1-ppos)
		if(ppos==0 & psdpos==0){csdpos <- 0}
		concpos <- pmax(-log(1-ppos)+c(-1,1)*qnorm(0.975)*sqrt(csdpos^2*(1+1/length(NTCxN))+(1/(1-ppos)-1)/ndrop),0)
	}
	else{ pirci <- c(0,1)
		conci <- c(0,Inf)
		cipos <- c(0,ndrop)
		cineg <- c(0,ndrop)
		cirain <- c(0,ndrop)
		psd <- NA
		csd <- NA
		psdpos <- NA
		csdpos <- NA
		concpos <- c(0,Inf)}

	try(conc <- rbind(c(1-px,-log(px)*1000/vol,ppos,-log(1-ppos)*1000/vol),
				c(psd,csd*1000/vol,psdpos,csdpos*1000/vol),
				c(1-pirci[2],conci[1]*1000/vol,cipos[1]/ndrop,concpos[1]*1000/vol),
				c(1-pirci[1],conci[2]*1000/vol,cipos[2]/ndrop,concpos[2]*1000/vol)))
	try(colnames(conc) <- c("prob_robust","conc_robust","prob_thres","conc_thres"))
	try(rownames(conc) <- c("est","sd","CI LB","CI UB"))
		# column 1: probability to be positive with robust estimator
		# column 2: concentration with robust estimator
		# column 3: probability to be positive with 1% threshold
		# column 4: concentration with 1% threshold
		# row 1: estimators
		# row 2: standard deviations
		# row 3: lower bound confidence intervals (including Poisson)
		# row 4: upper bound confidence intervals (including Poisson)

	try(tresh <- rbind(c(npos/ndrop,cipos/ndrop,npos,cipos),
				c(nneg/ndrop,cineg/ndrop,nneg,cineg),
				c(nrain/ndrop,cirain/ndrop,nrain,cirain)))
	try(colnames(tresh) <- c("prob","CI low p","CI high p","count","CI low n","CI high n"))
	try(rownames(tresh) <- c("pos","neg","rain"))

# convergence message
	time2 <- Sys.time()
	timerun <- time2-time1

#cat("Timetest\n")
#cat("Objects:",Ttime1-time1,"\n")
#cat("locfunc:",Ttime2-Ttime1,"\n")
#cat("isotonic:",Ttime3-Ttime2,"\n")
#cat("write output:",time2-Ttime3,"\n")

	cat("Estimation",modout$name,"completed in",timerun,attr(timerun,"units"),". \n")

# Print the results

cat("\n")
cat("Well",modout$name,":\n")
cat("\n")
cat(ndrop,"partitions\n")
cat("Robust estimator:\n")
cat("  ",round(ndrop*conc[1,1]),"positive, p =",round(conc[1,1],4),"( CI: [",round(conc[3,1],4),",",round(conc[4,1],4),"] )\n")
cat("  ",round(conc[1,2],1),"copies/ml ( CI: [",round(conc[3,2],1),",",round(conc[4,2],1),"] )\n")
cat("Threshold estimator:\n")
cat("  ",round(ndrop*conc[1,3]),"positive, p =",round(conc[1,3],4),"( CI: [",round(conc[3,3],4),",",round(conc[4,3],4),"] )\n")
cat("  ",round(conc[1,4],1),"copies/ml ( CI: [",round(conc[3,4],1),",",round(conc[4,4],1),"] )\n")
cat("CI's for concentration include Poisson variability.\n")
cat("__________________________________________________\n")


rlist <- list(conc=conc,tresh=tresh,fits=prl,data=data,droppi=probx,dropci=cbind(probciL,probciU),reppi=probr,densfits=datafit,ntcfits=ntcfit,pifits=pifits)
return(rlist)
}



#######################
######   Shell   ######
#######################

if(any(duplicated(names(datalist)))){
stop("Duplicated names in data, please provide unique names for each reaction")}

### Step 1: Identify negative cluster centre per well

modout <- lapply(datalist,modeRet)
names(modout) <- names(datalist)
for(i in 1:length(datalist)){modout[[i]]$name <- names(datalist)[[i]]}

NTCbin <- rep(0,length(datalist))
if(!is.null(NTC)){
	for(i in 1:length(NTC)){
		if(any(names(modout)==NTC[i])){
			NTCbin[which(names(datalist)==NTC[i])] <- 1
			modout[[which(names(datalist)==NTC[i])]]$mode <- modeRet(datalist[[which(names(datalist)==NTC[i])]],NTC=T)
}}}
cat("Step 1: Identify negative cluster centre per well - completed\n")

### Step 2: Check cluster centres with plate info

allmodes <- sapply(modout,function(x){x$mode})
platmod <- c(median(allmodes,na.rm=T),huber(allmodes[!is.na(allmodes)])$s)
if(length(datalist) > 2){modplat <- lapply(modout,checkmod1)
cat("Step 2: Check centres with plate info - completed\n")}
else{modplat <- modout
cat("Step 2: Check centres with plate info - skipped\n")}

### Step 3: Fit NTC's

	# create univariate NTC objects, all NTC's in one object
NTCx <- NULL
NTCxN <- NULL
if(!is.null(NTC)){
	NTCx <- list()
	for(i in 1:length(NTC)){
		if(any(names(datalist)==NTC[i])){
			NTCx <- append(NTCx,list(datalist[[which(names(datalist)==NTC[i])]]-modplat[[which(names(modplat)==NTC[i])]]$mode))
				# NTC fluorescence minus null mode
	}}
		# if there are NTC's, the NTCx object includes all NTC channel data
	if(!is.null(NTCx)){
		NTCxN <- lapply(NTCx,NTCfunc)
			# fit the NTC null distribution
		cat("Step 3: NTC identification - completed\n")}
	else{cat("      Warning: no wells found with given NTC names.\n")
		cat("Step 3: NTC identification - skipped\n")}}
else{cat("Step 3: NTC identification - skipped\n")}

	# This is the real procedure, all calculations happen internally in this step
resultlist <- lapply(as.list(1:length(datalist)),function(i){Shellfunc(datalist[[i]],modplat[[i]])})
names(resultlist) <- names(datalist)
for(i in 1:length(datalist)){resultlist[[i]]$name <- names(datalist)[[i]]}
cat("Step 4: Estimation - complete\n")


return(resultlist)

}
