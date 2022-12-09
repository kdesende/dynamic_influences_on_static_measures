# v-ratio example script
#
# This script provides an example as how to estimate v-ratio based on an experiment with 1 difficulty level, 2 choice option, and N confidence levels
# For more information about v-ratio, read the paper: https://pubmed.ncbi.nlm.nih.gov/35864100/
#
# This code was written in R-4.1.2, using the following package versions:
# Rcpp v1.0.7
# RcppZiggurat v0.1.6
# DEoptim v2.2.6
# reshape v0.8.8
#
# Code written by Kobe Desender and Luc Vermeylen, 2/12/2022

# House Keeping -----------------------------------------------------------
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory to current
rm(list=ls())

# Required Libraries
library(Rcpp) # to source, compile and run C++ functions
library(RcppZiggurat) #required for additional random number generation functions
# Note that Rcpp requires a c++ compiler, so you might have to install this; if so have a look here: https://teuder.github.io/rcpp4everyone_en/020_install.html
library(DEoptim) # optimization algorithm
library(reshape)
sourceCpp("DDM_with_confidence_slow_fullconfRT.cpp") # this will give R access to the function that is used to simulate from the v-ratio model (written in rcpp to increase speed)

# Below is the  objective function that is used to fit the model. It requires three inputs:
# Params (in the following order):  
### v=drift rate (usual range 0-4)
### a=boundary (usual range .5-4)
### ter=non-decision time (usual range 0-2)
### vratio=the measure of interest; note that the scripts _directly_ estimates v-ratio (usual range 0-2)
### addmean=parameter to map confidence predictions to the empirical data
### addsd=parameter to map confidence predictions to the empirical data

# observations: the data you want to fit. Note that the following columns are required:
### "cor": accuracy, with 1=correct, 0=error
### "rt": reaction time, in seconds
### "rtcj": reaction time of the confidence judgments, in seconds,
### "cj": confidence judgment

#conf_levels: the number of confidence levels in your experiment. Defaults to 6, implying the levels 1 (=lowest confidence),2,3,4,5,6 (=highest confidence)

# returnFit: a logical indicating whether you want to return the fit (1) or predictions given these parameters (0)

chi_square_optim <- function(params, observations, conf_levels=6, returnFit){
  sourceCpp("DDM_with_confidence_slow_fullconfRT.cpp") # for parallel use, call the simulation script again in the function
  
  #Some defaults
  z <- .5 #starting point; note that the mode is accuracy coded, so doesn't allow for starting point estimation currently
  ntrials <- 5000 #number of simulated data points
  dt <- .001 #precision of the simulation
  sigma <- 1 #within-trial noise, fixed at 1 (note that this is sometimes fixed at .1 in the Ratcliff tradition)
  
  # 1. generate predictions from the given parameters
  names(params) <- c('v','a','ter','vratio','add_mean','add_sd')
  
  #we're using the empirically observed confidence RTs for the simulations, stitching them X times to reach sufficient simulations
  rtcj_simuls <- rep(observations$rtcj,ceiling(ntrials/length(observations$rtcj))) 
  #Simulate data 
  predictions <- data.frame(DDM_with_confidence_slow_fullconfRT(v=params['v'],a=params['a'],ter=params['ter'],z=z,ntrials=length(rtcj_simuls),s=sigma,dt=dt,t2distribution=rtcj_simuls,postdriftmod=params['vratio']))
  
  #rt=predicted RT, resp=predicted resp, cor=predicted accuracy, raw_evidence2 = evidence at the time when confidence is quantified, rt2=predicted rtcj, cj=predicted cj
  names(predictions) <- c('rt','resp','cor','raw_evidence2','rtcj','cj') 

  # 2. Linear scaling of confidence, to transform them from a.u. onto a 1-6 scale (note; scale can be changed by varying conf_levels)
  predictions$cj <- (predictions$cj / params['add_sd']) + params['add_mean'] 
  predictions$cj[predictions$cj<1] <- 1
  predictions$cj[predictions$cj>conf_levels] <- conf_levels
  
  #only round numbers
  predictions$cj <- round(predictions$cj)
  
  # if we're only simulating data, return the predictions
  if(returnFit==0){ 
    return(predictions[,c('rt','cor','cj','rtcj')])
    
  # If we are fitting the model, now we'll compare predictions to the observations and compute the fit between both 
  }else{ 
    
    # First, separate the data in correct and error trials
    c_observed <- observations[observations$cor == 1,]
    e_observed <- observations[observations$cor == 0,]
    
    # same for the predictions
    c_predicted <- predictions[predictions$cor == 1,]
    e_predicted <- predictions[predictions$cor == 0,]

    # Now, get the quantile RTs on the "observed data" for correct and error distributions separately (for quantiles .1, .3, .5, .7, .9)
    c_quantiles <- quantile(c_observed$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    e_quantiles <- quantile(e_observed$rt, probs = c(.1,.3,.5,.7,.9), names = FALSE)
    
    # to combine correct and incorrect we scale the expected interquantile probability by the proportion of correct and incorect respectively
    prop_obs_c <- dim(c_observed)[1] / dim(observations)[1]
    prop_obs_e <- dim(e_observed)[1] / dim(observations)[1]
    
    c_obs_proportion = prop_obs_c * c(.1, .2, .2, .2, .2, .1)
    e_obs_proportion = prop_obs_e * c(.1, .2, .2, .2, .2, .1)
    obs_props <- c(c_obs_proportion,e_obs_proportion)
    
    # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data (scale by N?)
    c_pred_proportion <- c(
      sum(c_predicted$rt <= c_quantiles[1]),
      sum(c_predicted$rt <= c_quantiles[2]) - sum(c_predicted$rt <= c_quantiles[1]),
      sum(c_predicted$rt <= c_quantiles[3]) - sum(c_predicted$rt <= c_quantiles[2]),
      sum(c_predicted$rt <= c_quantiles[4]) - sum(c_predicted$rt <= c_quantiles[3]),
      sum(c_predicted$rt <= c_quantiles[5]) - sum(c_predicted$rt <= c_quantiles[4]),
      sum(c_predicted$rt > c_quantiles[5])
    ) / dim(predictions)[1]
    
    e_pred_proportion <- c(
      sum(e_predicted$rt <= e_quantiles[1]),
      sum(e_predicted$rt <= e_quantiles[2]) - sum(e_predicted$rt <= e_quantiles[1]),
      sum(e_predicted$rt <= e_quantiles[3]) - sum(e_predicted$rt <= e_quantiles[2]),
      sum(e_predicted$rt <= e_quantiles[4]) - sum(e_predicted$rt <= e_quantiles[3]),
      sum(e_predicted$rt <= e_quantiles[5]) - sum(e_predicted$rt <= e_quantiles[4]),
      sum(e_predicted$rt > e_quantiles[5])
    ) / dim(predictions)[1]
    pred_props <- c(c_pred_proportion,e_pred_proportion)
    
    # avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
    pred_props[pred_props==0] <- .0000001
    
    # Now, do the same for confidence
    obs_props_cj <- rep(NA,conf_levels*2)
    for(i in 1:conf_levels){
      obs_props_cj[i] <- sum(c_observed$cj==i)/dim(observations)[1]
      obs_props_cj[i+conf_levels] <- sum(e_observed$cj==i)/dim(observations)[1]
    }

    # now, get the proportion of responses that fall between the observed quantiles when applied to the predicted data (scale by N?)
    pred_props_cj <- rep(NA,conf_levels*2)
    for(i in 1:conf_levels){
      pred_props_cj[i] <- sum(c_predicted$cj==i)/dim(predictions)[1]
      pred_props_cj[i+conf_levels] <- sum(e_predicted$cj==i)/dim(predictions)[1]
    }

    # avoid zeros in the the data (because of division by predictions for chi square statistic) -> set to small number
    pred_props_cj[pred_props_cj==0] <- .0000001
    
    # Combine the quantiles for rts and cj
    obs_props <- c(obs_props,obs_props_cj)
    pred_props <- c(pred_props,pred_props_cj)
    
    # calculate chi square
    chiSquare = sum( ( (obs_props - pred_props) ^ 2) / pred_props )
    
    if (is.na(chiSquare)) {
      chiSquare = 1000000 # bad but quick solution to the NA's in objective function error... probably occurs when there are no errors...
    }
    return(chiSquare)
  }
}

# 1. The chi_square_optim function can be used for simulating data under the v-ratio model (returnFit=0)
# or it can be used to estimate v-ratio based on empirical data
# First, we will simulate some data under the v-ratio model, assuming some reasonable parameters
real_params <- c(1,1,.2,.75,4,1)
names(real_params) <- c("v","a","ter","vratio","add_mean","add_sd")

print(paste('Let us simulate data assuming the following parameters:',names(real_params),'=',real_params))

# In order to simulate data (under this version) of the v-ratio model, we need to know empirical confidence RTs
# For practical purpose, we just simulate some trials from a lognorm distribution and assume these reflect confidence RTs  
observations <- data.frame(matrix(NA,nrow=890));observations$rtcj <- rlnorm(890, log(1.5), log(1.5))

#how many confidence levels do you want in your simulations (note; 6 implies a 1-2-3-4-5-6 scale)?
conf_levels = 6

#Simulate data under the v-ratio model with these parameters
df = chi_square_optim(real_params,observations,conf_levels,returnFit=0) 

#some quick checks to reassure the data are meaningful
print(paste('average accuracy=',mean(df$cor)))
hist(df$rt[df$cor==1],col=rgb(0,1,0,.5),breaks=50,main="histogram of reaction times",xlab="RT")
hist(df$rt[df$cor==0],col=rgb(1,0,0,.5),add=T,breaks=50)
legend('bottom',legend=c("Corrects","Errors"),col=c("green","red"),box.lty=0,lty=1,inset=.1)
table(df$cj,df$cor) #distribution of confidence judgments, seperately for corrects and errors 

# 2. Now let's fit these data (assuming they are real data) to find the underlying parameters
# Usually, we don't know the true underlying parameters (such as when fitting real data),
# but in this case we know the ground truth, so we can check how good the model recovers the true parameters

# some fitting details 
itermax <- 1000 #ideally high enough (~1000, depending on convergence)
trace <- 50 #show the fitting output each X iterations

# parameter bounds for each of the parameters
v_range <- c(0,4)
a_range <- c(.5,4)
ter_range <- c(.1,1)
v_ratio_range <- c(0,2.5)
add_mean_range <- c(0,10) #note, these upper bounds depend a lot on the confidence scale used
add_sd_range <- c(0,10) #note, these upper bounds depend a lot on the confidence scale used

# fit the model to the simulated data
fit <- DEoptim(chi_square_optim, # function to optimize
               lower=c(v_range[1],a_range[1],ter_range[1],v_ratio_range[1],add_mean_range[1],add_sd_range[1]), # v,a,ter,vratio,add_mean,add_sd
               upper=c(v_range[2],a_range[2],ter_range[2],v_ratio_range[2],add_mean_range[2],add_sd_range[2]),
               observations=df,
               conf_levels=conf_levels,
               returnFit=1,
               control=list(itermax=itermax, trace=trace, parallelType=1, packages=c("Rcpp")))
fitted_params <- round(fit$optim$bestmem, 4)
names(fitted_params) <- c("v","a","ter","vratio","add_mean","add_sd")

#Voila, you have your parameters :)
print(paste('for this participant, v=',fitted_params['v'],'a=',fitted_params['a'],'ter=',fitted_params['ter'],', and most importantly, v-ratio =',fitted_params['vratio']))

#Compare these to the generative parameters
print(paste('the true parameters were the following: v=',real_params['v'],'a=',real_params['a'],'ter=',real_params['ter'],', andv-ratio =',real_params['vratio']))

#Finally, generate some simulations based on these parameters, and plot these on top of the empirical data to assess the fit
Simuls  = chi_square_optim(fitted_params,observations=df,conf_levels,returnFit=0)

#plotting empirical data and model fit on top of it
par(mfrow=c(1,2),mar=c(5.1,5.1,4.1,2.1))
tempC <- hist(df$rt[df$cor==1],breaks=seq(0,6.2,.1),xlim=c(0,2),prob=F,col=rgb(0,1,0,.25),border="white",ylab="Frequency",xlab="Reaction times (s)",cex.lab=2, cex.main=1.5, cex.axis=1.5,main="")
tempE <- hist(df$rt[df$cor==0],breaks=seq(0,6.2,.1),prob=F,add=T,col=rgb(1,0,0,.25),border='white')
Cors <- hist(Simuls$rt[Simuls$cor==1],breaks=seq(0,20,.1),plot=F)
Errs <- hist(abs(Simuls$rt[Simuls$cor==0]),breaks=seq(0,20,.1),plot=F)
lines(Cors$counts/(sum(Cors$counts)/sum(tempC$counts))~Cors$mids,type='l',col='green',lwd=3)
lines(Errs$counts/(sum(Errs$counts)/sum(tempE$counts))~Errs$mids,type='l',col='red',lwd=3)
#same for confidence
df$cj <- factor(df$cj,levels=1:conf_levels);Simuls$cj <- factor(Simuls$cj,levels=1:conf_levels)
df$ones <- 1; tempCJ <- with(df,aggregate(ones,by=list(cj=cj,cor=cor),sum));tempCJ <- cast(tempCJ,cor~cj);tempCJ[is.na(tempCJ)] <- 0
Simuls$ones <- 1; tempCJsim <- with(Simuls,aggregate(ones,by=list(cj=cj,cor=cor),sum));tempCJsim <- cast(tempCJsim,cor~cj);tempCJsim[is.na(tempCJsim)] <- 0
tempC <- barplot(as.matrix(tempCJ[tempCJ$cor==1,2:(conf_levels+1)]),ylim=c(0,max(max(tempCJ),max(tempCJsim))),col=rgb(0,1,0,.25),border='white',ylab="Frequency",xlab="Confidence",cex.lab=2, cex.main=1.5, cex.axis=1.5,main="")
tempE <- barplot(as.matrix(tempCJ[tempCJ$cor==0,2:(conf_levels+1)]),add=T,col=rgb(1,0,0,.25),border='white',axes=F)
lines(tempC,tempCJsim[tempCJsim$cor==1,2:(conf_levels+1)]/(dim(Simuls)[1]/dim(df)[1]),type='p',col='green',pch=4,lwd=3)
lines(tempC,tempCJsim[tempCJsim$cor==0,2:(conf_levels+1)]/(dim(Simuls)[1]/dim(df)[1]),type='p',col='red',pch=4,lwd=3)



