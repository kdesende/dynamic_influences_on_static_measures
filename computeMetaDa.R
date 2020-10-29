##########################################
### Code to compute meta-da in R #########
##########################################

# computeMetaDa computes meta-da based on the method described by Maniscalco, B., & Lau, H. (2012).
# A signal detection theoretic approach for estimating metacognitive sensitivity from confidence ratings.
# Consciousness and Cognition, 21, 420-430. doi: 10.1016/j.concog.2011.09.021

# Arguments:
#   ratings: a factor with levels corresponding to the rating categories,
#            ordered from low to high
#   stimulus: a factor, levels corresponding to stimulus identities
#   correct: a vector with 1 indicating correct responses and 0 incorrect responses
#   distr: What distributions of the evidence should be assumed. Default "norm", uses the normal distribution,
#          "logis" assumes the logistic distribution, and "gumbel" the distribution of smallest extremes
#   addConstant:  a logical value indicating if a constant of .5  should be added to each cell
#   nIterations: maximum number of iterations before the fitting procedure should be aborted
#
# Returns: a list with the elements
#   Da:  estimated sensitivity from objective discrimination responses
#   metaDa: estimated sensitivity from rating data
#   SdRatio: Ratio of the standard deviations of the signal created by the two stimuli
#   criteria: estimated criteria with respect to the optimal observer signal detection model
#   logLikelihood: logLikelihood of the best fit

# Our code avoids the slow linear inequality constraints of the original method, but fits the logs of parameters instead


computeMetaDa   <- function(ratings, stimulus, correct,
                              distr= "norm", addConstant = TRUE, nIterations = 10^3) {
  if(!is.factor(ratings)) stop ("ratings should be a factor!")
  if(!is.factor(stimulus )|| length(levels(stimulus)) != 2) {
    stop("stimulus should be a factor with 2 levels")
  }
  if(!all(correct %in% c(0,1))) stop("correct should be 1 or 0")

  pfun <- switch(distr, norm = pnorm, logis = plogis,
                 gumbel = function(x, location, scale) exp(-exp((location-x)/scale)))
  qfun <- switch(distr, norm = qnorm, logis = qlogis,
                 gumbel = function(x) -log(-log(x)))

  nRatings <-  length(levels(ratings))
  nCriteria <- nRatings * 2 - 1
  abs_corrects <-  table(ratings[correct == 1], stimulus[correct == 1])
  abs_errors   <-  table(ratings[correct == 0], stimulus[correct == 0])

  if (addConstant){
    abs_corrects <- abs_corrects + .5
    abs_errors   <- abs_errors + .5
  }

  abs_S1 <- c(rev(abs_errors[,2]),abs_corrects[,2])
  ratingHrs <- qfun(1 - cumsum(abs_S1)/sum(abs_S1))
  abs_S2 <-  c(rev(abs_corrects[,1]), abs_errors[,1] )
  ratingFrs <-  qfun(1 - cumsum(abs_S2)/sum(abs_S2))
  finits <- is.finite(ratingHrs) & is.finite(ratingFrs)
  ratingHrs <- as.vector(ratingHrs[finits])
  ratingFrs <- as.vector(ratingFrs[finits])

  nC_rS1 <- rev(as.vector(abs_corrects[,1]))
  nI_rS1 <- rev(as.vector(abs_errors[,2]))
  nC_rS2 <- as.vector(abs_corrects[,2])
  nI_rS2 <- as.vector(abs_errors[,1])

  meta_d1 <- ratingHrs[nRatings] - ratingFrs[nRatings]
  cs_1 <- -.5 * ( ratingHrs  + ratingFrs)
  initials <- c(meta_d1, 0,  cs_1[1], log(diff(cs_1)))

  optimLL <- function(parameters, nC_rS1,nI_rS1, nC_rS2,nI_rS2,nRatings, pfun){
    S1mu <- -parameters[1]/2
    S2mu <- parameters[1]/2
    S1sd <- 1
    S2sd <- 1/exp(parameters[2])
    t2c1x <-  c(-Inf, cumsum(c(parameters[3], exp(parameters[4:length(parameters)]))), Inf)
    t1c <- t2c1x[nRatings+1]
    prC_rS1 <- (pfun(t2c1x[2:(nRatings+1)],S1mu,S1sd) -
                pfun(t2c1x[1:nRatings],S1mu,S1sd) ) /pfun(t1c,S1mu,S1sd)
    prI_rS1 <- (pfun(t2c1x[2:(nRatings+1)],S2mu,S2sd) -
                pfun(t2c1x[1:nRatings],S2mu,S2sd) ) / pfun(t1c,S2mu,S2sd)
    prC_rS2 <- ((1- pfun(t2c1x[(nRatings+1):(nRatings*2)],S2mu,S2sd)) -
               (1- pfun(t2c1x[(nRatings+2):(nRatings*2+1)],S2mu,S2sd)) ) /
               ( 1 - pfun(t1c,S2mu,S2sd))
    prI_rS2 <- ((1- pfun(t2c1x[(nRatings+1):(nRatings*2)],S1mu,S1sd)) -
               (1- pfun(t2c1x[(nRatings+2):(nRatings*2+1)],S1mu,S1sd)) ) /
                (1 - pfun(t1c,S1mu,S1sd))
    logL <-  -sum(nC_rS1*log(prC_rS1),nI_rS1*log(prI_rS1),nC_rS2*log(prC_rS2),nI_rS2*log(prI_rS2))
    logL
  }

  fit <- optim(par = initials, f = optimLL, gr = NULL,
                nC_rS1 = nC_rS1, nI_rS1 = nI_rS1, nC_rS2 = nC_rS2, nI_rS2 = nI_rS2, nRatings = nRatings,pfun,
                control = list(maxit = nIterations))
  s <- exp(fit$par[2])

 criteria <- cumsum(c(fit$par[3], exp(fit$par[4:length(fit$par)])))
 result <- list(Da = meta_d1 * s * sqrt(2/(1 + s^2)),
                metaDa = fit$par[1] * s * sqrt(2/(1 + s^2)),
                SdRatio = s,
                Criteria = criteria,
                logLikelihood = -fit$value
                )

  result
}

### Examples

stimulus <- factor(rep(c("A", "B"), 10))
ratings <- factor(c(rep(1:3,4), rep(3,times=8)))
correct <- c(rep(0,4), rep(1,6), rep(0,2), rep(1,8))

computeMetaDa(stimulus = stimulus, ratings = ratings, correct = correct)
computeMetaDa(stimulus = stimulus, ratings = ratings, correct = correct, distr = "logis")
computeMetaDa(stimulus = stimulus, ratings = ratings, correct = correct, distr = "gumbel")


computeMetaDa(stimulus = stimulus, ratings = ratings, correct = correct,
             addConstant = FALSE)


