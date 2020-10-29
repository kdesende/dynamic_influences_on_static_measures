################
#1. Simulations 
rm(list=ls())
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
# to source, compile and run C++ functions
library(Rcpp) 
sourceCpp("DDM_with_confidence_slow.cpp") # this will give R access to the DDM_with_confidence_slow function 
#Load code to compute metad
source('computeMetaDa.R')
#Load the heatmap reprenseint p(correct|evt,t,x)
pcor <- as.matrix(read.table('pcor_5secs.txt',header=T)) #sigma=.1
pcorVector <- as.vector(pcor)
evUp <- seq(0,1,by=.01);


#1. Simulate data with random selection of parameters
N<-100 #N participants
samples=500 #ntrials,multiplied by 2 below
t2 <- 1 #post-decision time
z <- .5 #starting point
sigma = 1 #within-trial noise
for(sub in 1:N){
  drift <- runif(1,0,2.5);drift <- c(drift,-drift) #to estimate bias, we need a positive and negative drift
  ter <- runif(1,.2,.6) #variation in non decision time
  bound = runif(1,.5,3)  #variation in boundary
  vratio <- runif(1,0,1.5) #variation in vratio (i.e., individual differences in metacognition)
  
  out1 <- data.frame(DDM_with_confidence_slow(v=drift[1],a=bound,ter=ter,z=z,ntrials=samples,s=sigma,dt=.001,t2time=t2,postdriftmod=vratio))
  out1$drift <- drift[1]
  out2 <- data.frame(DDM_with_confidence_slow(v=drift[2],a=bound,ter=ter,z=z,ntrials=samples,s=sigma,dt=.001,t2time=t2,postdriftmod=vratio))
  out2$drift <- drift[2]
  
  out <- rbind(out1,out2);names(out) <- c('rt','resp','cor','evidence2','rt2','cj','drift')
  
  #Convert cj into p(correct), so first flip
  out$evidence2[out$resp==1] <- out$evidence2[out$resp==1] - .5*bound #subtract the starting point
  out$evidence2[out$resp==-1] <- (-out$evidence2[out$resp==-1]) + .5*bound
  
  #match to the heatmap
  out$evidence2 <- out$evidence2*.1 #the heatmap is created using sigma=.1
  out$closest_evdnc2 <- match.closest(abs(out$evidence2),evUp)
  out$outrt2 <- out$rt2;
  out$outrt2[out$outrt2>5] <- 5 #heatmap doesn't go higher
  out$outrt2 <- out$outrt2*(dim(pcor)[1]/5) #scale with the heatmap, between 0 and 2000
  out$conft2 <- pcorVector[(out$closest_evdnc2-1)*dim(pcor)[1]+round(out$outrt2)]
  out$conft2[out$resp==1&out$evidence2<0] <- out$conft2[out$resp==1&out$evidence2<0] - 2*(out$conft2[out$resp==1&out$evidence2<0]-.5)
  out$conft2[out$resp==-1&out$evidence2<0] <- out$conft2[out$resp==-1&out$evidence2<0] - 2*(out$conft2[out$resp==-1&out$evidence2<0]-.5)
  
  #Use heatmap confidence instead of raw evidence
  out$cj <- out$conft2
  
  out$bound <- bound;out$starting_point <- z;out$within_trial_noise <- sigma;out$vratio <- vratio;out$rtconf <- t2;out$ter <- ter
  
  out$cj_quant <- as.numeric(cut(out$cj,5))
  out$sub <- sub
  
  #save to df
  if(sub==1){
    Data <- out
  }else{
    Data <- rbind(Data,out)
  }
}

#Show some sanity checks
par(mfrow=c(1,2))
cjAcc <- with(Data,aggregate(cor,by=list(confidence=cj_quant,sub=sub),mean));cjAcc<-cast(cjAcc,sub~confidence)
plot(colMeans(cjAcc,na.rm=T),type='b',xlab='confidence (divided into 6 bins)',ylab='accuracy',frame=F,ylim=c(0,1));error.bar(1:5,colMeans(cjAcc,na.rm=T),colSds(as.matrix(cjAcc),na.rm=T)/sqrt(N),length=0)
cjRT <- with(Data,aggregate(rt,by=list(confidence=cj_quant,sub=sub),median));cjRT<-cast(cjRT,sub~confidence) #the model doesn't predict this with a single level of coherence
plot(colMeans(cjRT,na.rm=T),type='b',xlab='confidence (divided into 6 bins)',ylab='median RT',frame=F,ylim=c(.5,1.2));error.bar(1:5,colMeans(cjRT,na.rm=T),colSds(as.matrix(cjRT),na.rm=T)/sqrt(N),length=0)

#Compute meta-d
d <- rep(NA,N);metad <- rep(NA,N);mratio<-rep(NA,N)
for(i in 1:N){
  temp <- subset(Data,sub==i)
  temp$cresp <- sign(temp$drift)
  temp$cj <- as.numeric(cut(temp$cj,4))
  fit <- computeMetaDa(as.factor(temp$cj),as.factor(temp$cresp),as.factor(temp$cor))
  d[i] <- fit$Da;metad[i] <- fit$metaDa;mratio[i]<-metad[i]/d[i]
}
x <- data.frame(cbind(d,metad,mratio));names(x) <- c("d","meta-d","mratio")
par(mfrow=c(1,1))
stripchart(x,at=1:3,jitter=.03, method = 'jitter',pch=21,vertical = TRUE,bg="grey",col='white',cex=1.5,ylim=range(x),frame=F)

Data$drift_abs <- abs(Data$drift)
drift <- with(Data,aggregate(drift_abs,by=list(sub),mean))$x
ter <- with(Data,aggregate(ter,by=list(sub),mean))$x
bound <- with(Data,aggregate(bound,by=list(sub),mean))$x
confidence_rt <- with(Data,aggregate(rtconf,by=list(sub),mean))$x
vratio <- with(Data,aggregate(vratio,by=list(sub),mean))$x

#Figure for the paper
jpeg(file="figure1_model_fits.jpg", res=600, width=400*3*(350/72), height=400*(350/72))
par(mfrow=c(1,3))
plot(mratio~vratio,frame=F,pch=21,cex=2,cex.lab=1.5,col="white",bg="black",ylab="M-ratio",xlab="vratio");cor.test(mratio,vratio);abline(lm(mratio~vratio),lty=2,lwd=2)
plot(mratio~bound,frame=F,pch=21,cex=2,cex.lab=1.5,col="white",bg="black",ylab="M-ratio",xlab="Decision boundary");cor.test(mratio,bound);abline(lm(mratio~bound),lty=2,lwd=2)
plot(vratio~bound,frame=F,pch=21,cex=2,cex.lab=1.5,col="white",bg="black",ylab="vratio",xlab="Decision boundary");cor.test(vratio,bound);abline(lm(vratio~bound),lty=2,lwd=2)
dev.off()  

#Table showing all values for the paper
library(Hmisc);library(grid);library(gridExtra)
x <- cbind(drift,ter,bound,vratio,mratio)
correlation_matrix<-rcorr(x, type="pearson")
R <- correlation_matrix$r # Matrix of correlation coeficients
p <- correlation_matrix$P # Matrix of p-value
mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
## build a new matrix that includes the correlations with their apropriate stars
Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
diag(Rnew) <- paste(diag(R), " ", sep="")
rownames(Rnew) <-   names(x) <- c('Drift rate','Non-decision time','Decision boundary','vratio','M-ratio')
colnames(Rnew) <-   names(x) <- c('Drift rate','Non-decision time','Decision boundary','vratio','M-ratio')
Rnew <- as.matrix(Rnew)
Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
Rnew <- as.data.frame(Rnew)
## remove last column and return the correlation matrix
Rnew <- cbind(Rnew[1:length(Rnew)-1])
jpeg(file="Manuscript dynamic metad/figure1_rawCors.jpg", res=600, width=400*3*(350/72), height=400*(350/72))
grid.newpage();grid.table(Rnew)
dev.off()
