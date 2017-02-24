PlotPartClass2d <- function(PartClass2dobject,platename,res=120){
	
	pal0 <- colorRampPalette(c("blue","red"))(76)
	pal1 <- colorRampPalette(c("black","gray"))(76)


	plotfunc <- function(wellres){
	name <- wellres$name
	ntcs <- length(wellres$ntcfits$x)
#	colx <- pal0[sqrt(wellres$densfits$x[[1]])/max(sqrt(wellres$densfits$x[[1]]))*75+1]
#	coly <- pal0[sqrt(wellres$densfits$y[[1]])/max(sqrt(wellres$densfits$y[[1]]))*75+1]
	colx <- pal0[wellres$droppi[,1]*75+1]
	coly <- pal0[wellres$droppi[,2]*75+1]
	cols <- pal1[sqrt(pmax(wellres$droppi[,1]*wellres$droppi[,2],wellres$droppi[,1]*(1-wellres$droppi[,2]),(1-wellres$droppi[,1])*wellres$droppi[,2],(1-wellres$droppi[,1])*(1-wellres$droppi[,2])))*150+1-75]

	# plot 1 - partition probabilities x & y + info 1D
	filename1 <- paste(platename,"_",name,"_","full_results_1D",".png",sep="")
	png(filename1,width=3600*res/120,height=1800*res/120,res=res)
	layout(matrix(c(2,1,4,3),nrow=2),widths=c(1,1), heights=c(1,3))
	par(oma=rep(0,4))
	par(mar=c(5,5.5,0,1))
		# plot x, colours x
		plot(wellres$data[,1],wellres$droppi[,1],pch=".",main="",xlab="fluorescence intensity",ylab=expression(hat(p)[i0](x)),xaxs="i",ylim=c(0,1),xlim=c(xmin,xmax),cex.axis=res/60,cex.lab=res/60)
		try(points(wellres$data[,1],wellres$dropci$x[,1],pch="."))
		try(points(wellres$data[,1],wellres$dropci$x[,2],pch="."))
		try(points(wellres$data[,1],wellres$droppi[,1],col=colx,pch=20))
		text(textlocx,0.85,paste(round(wellres$conc$x[1,1],4),"negative"),cex=res/40)
		text(textlocx,0.8,paste("CI: [",round(wellres$conc$x[3,1],4),",",round(wellres$conc$x[4,1],4),"]"),cex=res/40)
		text(textlocx,0.75,paste(round(wellres$conc$x[1,2],1),"copies/ml"),cex=res/40)
		text(textlocx,0.7,paste("CI: [",round(wellres$conc$x[3,2],1),",",round(wellres$conc$x[4,2],1),"]"),cex=res/40)
		lines(c((xmin+xmax)/2,xmax-xdiff*0.03),c(0.65,0.65))
		text(textlocx,0.6,paste(wellres$tresh$x[1,4],"positive droplets ([",wellres$tresh$x[1,5],",",wellres$tresh$x[1,6],"])"),cex=res/60)
		text(textlocx,0.55,paste(wellres$tresh$x[2,4],"negative droplets ([",wellres$tresh$x[2,5],",",wellres$tresh$x[2,6],"])"),cex=res/60)
		text(textlocx,0.5,paste(wellres$tresh$x[3,4],"rain droplets ([",wellres$tresh$x[3,5],",",wellres$tresh$x[3,6],"])"),cex=res/60)
	par(mar=c(0,5.5,3,1))
		xhist <- hist(wellres$data[,1],breaks=seq(xmin,xmax,xdiff/102),plot=F)
		barplot(xhist$intensities,space=0,horiz=F,axes = FALSE,col="purple",xaxs="i",main=paste(platename,name,"results Channel 1"),cex.main=res/40)
	par(mar=c(5,5.5,0,1))
		# plot y, colours y
		plot(wellres$data[,2],wellres$droppi[,2],pch=".",main="",xlab="fluorescence intensity",ylab=expression(hat(p)[i0](x)),xaxs="i",ylim=c(0,1),xlim=c(ymin,ymax),cex.axis=res/60,cex.lab=res/60)
		try(points(wellres$data[,2],wellres$dropci$y[,1],pch="."))
		try(points(wellres$data[,2],wellres$dropci$y[,2],pch="."))
		try(points(wellres$data[,2],wellres$droppi[,2],col=coly,pch=20))
		text(textlocy,0.85,paste(round(wellres$conc$y[1,1],4),"negative"),cex=res/40)
		text(textlocy,0.8,paste("CI: [",round(wellres$conc$y[3,1],4),",",round(wellres$conc$y[4,1],4),"]"),cex=res/40)
		text(textlocy,0.75,paste(round(wellres$conc$y[1,2],1),"copies/ml"),cex=res/40)
		text(textlocy,0.7,paste("CI: [",round(wellres$conc$y[3,2],1),",",round(wellres$conc$y[4,2],1),"]"),cex=res/40)
		lines(c((ymin+ymax)/2,ymax-ydiff*0.03),c(0.65,0.65))
		text(textlocy,0.6,paste(wellres$tresh$y[1,4],"positive droplets ([",wellres$tresh$y[1,5],",",wellres$tresh$y[1,6],"])"),cex=res/60)
		text(textlocy,0.55,paste(wellres$tresh$y[2,4],"negative droplets ([",wellres$tresh$y[2,5],",",wellres$tresh$y[2,6],"])"),cex=res/60)
		text(textlocy,0.5,paste(wellres$tresh$y[3,4],"rain droplets ([",wellres$tresh$y[3,5],",",wellres$tresh$y[3,6],"])"),cex=res/60)
	par(mar=c(0,5.5,3,1))
		yhist <- hist(wellres$data[,2],breaks=seq(ymin,ymax,ydiff/102),plot=F)
		barplot(yhist$intensities,space=0,horiz=F,axes = FALSE,col="purple",xaxs="i",main=paste(platename,name,"results Channel 2"),cex.main=res/40)
	dev.off()

	# plot 2 - 2D information
	filename2 <- paste(platename,"_",name,"_","full_results_2D",".png",sep="")
	png(filename2,width=3000*res/120,height=3000*res/120,res=res)
	par(oma=rep(0,4))
	par(mfrow=c(2,2))
		# y prob colour
	par(mar=c(5.5,5.5,4.5,1))
		plot(wellres$data[,c(1,2)],pch=20,col=coly,main=paste(platename,name,"probabilities channel 2"),xlab="fluorescence intensity channel 1",ylab="fluorescence intensity channel 2",cex.axis=res/60,cex.lab=res/60,cex.main=res/40)
		
		# information
	par(mar=c(0,0,0,0))
		plot(10,10,xlim=c(0,1),ylim=c(0,1),axes = FALSE,bty="n",xlab="",ylab="")
		text(0.25,0.9,"Channel 1",cex=res/30)
		text(0.25,0.8,paste(round(wellres$conc$x[1,1],4),"negative"),cex=res/40)
		text(0.25,0.75,paste("CI: [",round(wellres$conc$x[3,1],4),",",round(wellres$conc$x[4,1],4),"]"),cex=res/40)
		text(0.25,0.65,paste(round(wellres$conc$x[1,2],1),"copies/ml"),cex=res/40)
		text(0.25,0.6,paste("CI: [",round(wellres$conc$x[3,2],1),",",round(wellres$conc$x[4,2],1),"]"),cex=res/40)
		lines(c(0.05,0.45),c(0.5,0.5))
		text(0.25,0.4,paste(wellres$tresh$x[1,4],"positive droplets ([",wellres$tresh$x[1,5],",",wellres$tresh$x[1,6],"])"),cex=res/60)
		text(0.25,0.35,paste(wellres$tresh$x[2,4],"negative droplets ([",wellres$tresh$x[2,5],",",wellres$tresh$x[2,6],"])"),cex=res/60)
		text(0.25,0.3,paste(wellres$tresh$x[3,4],"rain droplets ([",wellres$tresh$x[3,5],",",wellres$tresh$x[3,6],"])"),cex=res/60)
		text(0.75,0.9,"Channel 2",cex=res/30)
		text(0.75,0.8,paste(round(wellres$conc$y[1,1],4),"negative"),cex=res/40)
		text(0.75,0.75,paste("CI: [",round(wellres$conc$y[3,1],4),",",round(wellres$conc$y[4,1],4),"]"),cex=res/40)
		text(0.75,0.65,paste(round(wellres$conc$y[1,2],1),"copies/ml"),cex=res/40)
		text(0.75,0.6,paste("CI: [",round(wellres$conc$y[3,2],1),",",round(wellres$conc$y[4,2],1),"]"),cex=res/40)
		lines(c(0.55,0.95),c(0.5,0.5))
		text(0.75,0.4,paste(wellres$tresh$y[1,4],"positive droplets ([",wellres$tresh$y[1,5],",",wellres$tresh$y[1,6],"])"),cex=res/60)
		text(0.75,0.35,paste(wellres$tresh$y[2,4],"negative droplets ([",wellres$tresh$y[2,5],",",wellres$tresh$y[2,6],"])"),cex=res/60)
		text(0.75,0.3,paste(wellres$tresh$y[3,4],"rain droplets ([",wellres$tresh$y[3,5],",",wellres$tresh$y[3,6],"])"),cex=res/60)
				
		# original data
	par(mar=c(5.5,5.5,4.5,1))
		plot(wellres$rawdata,pch=20,main=paste(platename,name,"Raw data"),xlab="fluorescence intensity channel 1",ylab="fluorescence intensity channel 2",cex.axis=res/60,cex.lab=res/60,cex.main=res/40)

		# x prob colour
	par(mar=c(5.5,5.5,4.5,1))
		plot(wellres$data[,c(1,2)],pch=20,col=colx,main=paste(platename,name,"probabilities channel 1"),xlab="fluorescence intensity channel 1",ylab="fluorescence intensity channel 2",cex.axis=res/60,cex.lab=res/60,cex.main=res/40)
	dev.off()
	}


xmin <- min(sapply(PartClass2dobject,function(x){min(x$data[,1],na.rm=T)}),na.rm=T)
xmax <- max(sapply(PartClass2dobject,function(x){max(x$data[,1],na.rm=T)}),na.rm=T)
ymin <- min(sapply(PartClass2dobject,function(x){min(x$data[,2],na.rm=T)}),na.rm=T)
ymax <- max(sapply(PartClass2dobject,function(x){max(x$data[,2],na.rm=T)}),na.rm=T)
xdiff <- xmax-xmin
ydiff <- ymax-ymin
xmin <- xmin-xdiff/100
xmax <- xmax+xdiff/100
ymin <- ymin-ydiff/100
ymax <- ymax+ydiff/100
textlocx <- xmin+xdiff*0.75
textlocy <- ymin+ydiff*0.75

lapply(PartClass2dobject,plotfunc)

# plot all together - pifits
	filename0 <- paste(platename,"_","pi0fits",".png",sep="")
	nwell <- length(PartClass2dobject)
	png(filename0,width=1800*res/120,height=900*nwell*res/120,res=res)
	par(mfrow=c(nwell,2))
	par(mar=c(4.5,4.5,4.5,1))
	pifitting <- function(wellres){
		ntcs <- length(wellres$ntcfits$x)
		plot(c(-4,3),rep(log(1-wellres$conc$x[1,1]),2),xlim=c(-3,2),ylim=log(1-wellres$conc$x[1,1])+c(-0.5,0.5),ylab="log(pi0)",xlab="location (sigma's from mu)",cex.lab=res/60,cex.axis=res/60,main=paste(platename,wellres$name,"Channel 1"),cex.main=res/60)
	for(m in 1:ntcs){
		try(lines(seq(-3,2,0.1),log(wellres$pifits$x[[m]])))
		try(lines(seq(-1,0,0.1),log(wellres$pifits$x[[m]][21:31]),lwd=3))
		try(lines(seq(-1,0,0.1),rep(log(wellres$fits$pi[m,1]),11),lwd=1,col="red"))
		}
		ntcs <- length(wellres$ntcfits$y)
		plot(c(-4,3),rep(log(1-wellres$conc$y[1,1]),2),xlim=c(-3,2),ylim=log(1-wellres$conc$y[1,1])+c(-0.5,0.5),ylab="log(pi0)",xlab="location (sigma's from mu)",cex.lab=res/60,cex.axis=res/60,main=paste(platename,wellres$name,"Channel 2"),cex.main=res/60)
	for(m in 1:ntcs){
		try(lines(seq(-3,2,0.1),log(wellres$pifits$y[[m]])))
		try(lines(seq(-1,0,0.1),log(wellres$pifits$y[[m]][21:31]),lwd=3))
		try(lines(seq(-1,0,0.1),rep(log(wellres$fits$pi[m,2]),11),lwd=1,col="red"))
		}
	}
	lapply(PartClass2dobject,pifitting)
	dev.off()
}




