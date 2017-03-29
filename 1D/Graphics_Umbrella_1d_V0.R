PlotUmbrella1d <- function(Umbrella1dObject,platename,res=120){
	
	pal0 <- colorRampPalette(c("red","blue"))(76)

	plotfunc <- function(wellres){
	name <- wellres$name
	ntcs <- length(wellres$ntcfits)
	cols <- pal0[wellres$droppi*75+1]

	# plot - droplet probabilities + info
	filename0 <- paste(platename,"_",name,"_","full_results",".png",sep="")
	png(filename0,width=1800*res/120,height=1800*res/120,res=res)
	layout(matrix(c(2,1),nrow=2),widths=1, heights=c(1,3))
	#1
	par(mar=c(5,5.5,0,1))
		plot(sort(wellres$data),wellres$droppi[order(wellres$data)],type="l",main="",xlab="Fluorescence intensity",ylab=expression(paste(hat(p)[i0](x),"  (probability to be negative)")),xaxs="i",ylim=c(0,1),xlim=c(xmin,xmax),cex.axis=2,cex.lab=2)
		try(abline(h=0.8,col="dark gray",lwd=2))
		try(abline(h=0.05,col="dark gray",lwd=2))
		try(lines(sort(wellres$data),wellres$dropci[order(wellres$data),1],lty=2))
		try(lines(sort(wellres$data),wellres$dropci[order(wellres$data),2],lty=2))
		points(wellres$data,wellres$droppi,col=cols,pch=20)
		text(textloc,0.75,paste("Robust estimator:       "),cex=1.5)
		text(textloc,0.71,paste(round(wellres$conc[1,1],4),"positive"),cex=1.5)
		text(textloc,0.67,paste("CI: [",round(wellres$conc[3,1],4),",",round(wellres$conc[4,1],4),"]"),cex=1.5)
		text(textloc,0.63,bquote(.(round(wellres$conc[1,2],1))~copies / mu*l),cex=1.5)
		text(textloc,0.59,paste("CI: [",round(wellres$conc[3,2],1),",",round(wellres$conc[4,2],1),"]"),cex=1.5)
		text(textloc,0.51,paste("Threshold estimator:    "),cex=1.5)
		text(textloc,0.47,paste(round(wellres$conc[1,3],4),"positive"),cex=1.5)
		text(textloc,0.43,paste("CI: [",round(wellres$conc[3,3],4),",",round(wellres$conc[4,3],4),"]"),cex=1.5)
		text(textloc,0.39,bquote(.(round(wellres$conc[1,4],1))~copies / mu*l),cex=1.5)
		text(textloc,0.35,paste("CI: [",round(wellres$conc[3,4],1),",",round(wellres$conc[4,4],1),"]"),cex=1.5)
		text(textloc,0.92,paste(wellres$tresh[1,4],"positive partitions ([",wellres$tresh[1,5],",",wellres$tresh[1,6],"])"),cex=1.5)
		text(textloc,0.88,paste(wellres$tresh[2,4],"negative partitions ([",wellres$tresh[2,5],",",wellres$tresh[2,6],"])"),cex=1.5)
		text(textloc,0.84,paste(wellres$tresh[3,4],"rain partitions ([",wellres$tresh[3,5],",",wellres$tresh[3,6],"])"),cex=1.5)
	#2
	par(mar=c(0,5.5,3,1))
	xhist <- hist(wellres$data,breaks=seq(xmin,xmax,xdiff/102),plot=F)
	barplot(xhist$density,space=0,horiz=F,axes = FALSE,col="purple",xaxs="i",main=paste(platename,name,"results"),cex.main=3)
	dev.off()
	}

xmin <- min(sapply(Umbrella1dObject,function(x){min(x$data,na.rm=T)}),na.rm=T)
xmax <- max(sapply(Umbrella1dObject,function(x){max(x$data,na.rm=T)}),na.rm=T)
xdiff <- xmax-xmin
xmin <- xmin-xdiff/100
xmax <- xmax+xdiff/100
textloc <- xmin+xdiff*0.75

lapply(Umbrella1dObject,plotfunc)

# plot all together - pifits
	filename0 <- paste(platename,"_","pi0fits",".png",sep="")
	nwell <- length(Umbrella1dObject)
	png(filename0,width=900*res/120,height=900*nwell*res/120,res=res)
	par(mfrow=c(nwell,1))
	par(mar=c(4.5,4.5,4.5,1))
	pifitting <- function(wellres){
		ntcs <- length(wellres$ntcfits)
		plot(c(-4,3),rep(wellres$conc[1,1],2),xlim=c(-3,2),ylim=1-wellres$conc[1,1]+c(-0.2,0.2),ylab=expression(hat(p)[0]),xlab=expression(paste("location (",sigma," from ",mu,")")),cex.lab=2,cex.axis=2,main=paste(platename,wellres$name),cex.main=3,cex.axis=2,cex.lab=2)
	for(m in 1:ntcs){
		try(lines(seq(-3,2,0.1),wellres$pifits[[m]]))
		try(lines(seq(-1,0,0.1),wellres$pifits[[m]][21:31],lwd=3))
		try(lines(seq(-1,0,0.1),rep(wellres$fits$pi[m],11),lwd=1,col="red"))
		}
	}
	lapply(Umbrella1dObject,pifitting)
	dev.off()
}
