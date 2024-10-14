
# This script fits the 4-parameter ddm fits to the mystery data (2 stimulus
# conditions)

library("rtdists")
load("mystery.Rdata")

# Code for condition RTs is rt.[stimulus][response];
# re is correct/incorrect, not the response.
# Create a data frame, converting the RTs in ms to s
data <- data.frame(rt=rt,re=re,stim=stim)

# Pull out the different RTs for different conditions to evaluate
# model fit later, convert to seconds so we can use EZ-diffusion
# for the starting values.
rt.AA <- data$rt[data$re==1 & data$stim==0]/1000
rt.BB <- data$rt[data$re==1 & data$stim==1]/1000
rt.AB <- data$rt[data$re==0 & data$stim==0]/1000
rt.BA <- data$rt[data$re==0 & data$stim==1]/1000

# Starting values for RTs in seconds, from Wagenmakers' online
# EZ-diffusion parameter estimation
  # For A: a = .694, nu.A = .675,t0 = .414
  # For B: a = .733, nu.B = .936,t0 = .388

# I used the commented out starting values above in the first NM iteration, 
# but then I reran the code several times, using the recovered "best" parameter
# values as starting values for the next NM iteration.  I performed the NM
# optimization about 3 times (IIRC).
a <- .205 # log scale
z <- -.514 # will give z = a/2
nu.A <- 1.459 # log scale
nu.B <- .0193 # negative log scale
t0 <- -1.1137 # log scale

a = log(.710)
z <- 0
nu.A <- log(.675)
nu.B <- log(.936)
t0 <- log(.401)
start.pars <- c(a,z,nu.A,nu.B,t0)

# DDM likelihood computation
mle.diffusion <- function(pars,data) {
  # Reparameterize
  a <- exp(pars[1])
  z <- a/(1 + exp(-pars[2]))
  nu.A <- exp(pars[3])
  nu.B <- exp(pars[4])
  t0 <- exp(pars[5])
  
  # rt.[stimulus][response]
  # re is correct/incorrect, not the response
  rt.AA <- data$rt[data$re==1 & data$stim==0]/1000
  rt.BB <- data$rt[data$re==1 & data$stim==1]/1000
  rt.AB <- data$rt[data$re==0 & data$stim==0]/1000
  rt.BA <- data$rt[data$re==0 & data$stim==1]/1000
  
  # Compute the likelihoods for each condition.  The "A" response is
  # the upper boundary, "B" is the lower (0) boundary
  ll.AA <- ddiffusion(rt.AA,response="upper",a,nu.A,t0,z)
  ll.BB <- ddiffusion(rt.BB,response="upper",a,nu.B,t0,a-z)
  ll.AB <- ddiffusion(rt.AB,response="lower",a,nu.A,t0,z)
  ll.BA <- ddiffusion(rt.BA,response="lower",a,nu.B,t0,a-z)
  
  # Catch any 0s and replace them with something tiny
  ll.AA[ll.AA<=1e-200] <- 1e-200
  ll.BB[ll.BB<=1e-200] <- 1e-200
  ll.AB[ll.AB<=1e-200] <- 1e-200
  ll.BA[ll.BA<=1e-200] <- 1e-200
  
  # Compute the total log likelihood
  ll <- -(sum(log(ll.AA)) + sum(log(ll.BB)) + sum(log(ll.AB)) + sum(log(ll.BA)))
  return(ll)
  }

 # Maximize the likelihood
 fit.ddm <- optim(start.pars,mle.diffusion,data=data,control=list(maxit=10000))
 
 # Recover parameters
  a <- exp(fit.ddm$par[1])
  z <- a/(1 + exp(-fit.ddm$par[2]))
  nu.A <- exp(fit.ddm$par[3])
  nu.B <- exp(fit.ddm$par[4])
  t0 <- exp(fit.ddm$par[5])
 
 # Visualize the fits  
 par(mfrow=c(2,2))
 
 # Note that because ddiffusion returns the *joint* pdf and not the pdf conditioned on 
 # the response, we need to divide by the integral of the joint pdf so that the
 # predicted density drawn on the histograms integrates to 1.
 times <- seq(0,1.1,length=200)
 hist(rt.AA,breaks=20,freq=FALSE)
  lines(times,ddiffusion(times,response="upper",a,nu.A,t0,z)/
              pdiffusion(10*max(times),response="upper",a,nu.A,t0,z),col="purple4")  
 hist(rt.BB,breaks=20,freq=FALSE)
  lines(times,ddiffusion(times,response="upper",a,nu.B,t0,a-z)/
          pdiffusion(10*max(times),response="upper",a,nu.B,t0,a-z),col="purple4")
 hist(rt.AB,breaks=20,freq=FALSE)
  lines(times,ddiffusion(times,response="lower",a,nu.A,t0,z)/
          pdiffusion(10*max(times),response="lower",a,nu.A,t0,z),col="purple4")
 hist(rt.BA,breaks=20,freq=FALSE)
   lines(times,ddiffusion(times,response="lower",a,nu.B,t0,a-z)/
           pdiffusion(10*max(times),response="lower",a,nu.B,t0,a-z),col="purple4")

# If you want to rerun the optimization
   start.pars <- fit.ddm$par
   
# Maybe some QQ-plots, just for fun:
   ranks <- seq(.05,.95,by=.1)
   qs.AA <- quantile(rt.AA,ranks)
   qs.BB <- quantile(rt.BB,ranks)
   qs.AB <- quantile(rt.AB,ranks)
   qs.BA <- quantile(rt.BA,ranks)
   
   q.AA.model <- qdiffusion(ranks,response="upper",a,nu.A,t0,z,scale_p=TRUE)
   q.BB.model <- qdiffusion(ranks,response="upper",a,nu.B,t0,a-z,scale_p=TRUE)
   q.AB.model <- qdiffusion(ranks,response="lower",a,nu.A,t0,z,scale_p=TRUE)
   q.BA.model <- qdiffusion(ranks,response="lower",a,nu.B,t0,a-z,scale_p=TRUE)
   
   par(mfrow=c(2,2))
   plot(q.AA.model,qs.AA,xlab="Diffusion Quantiles",ylab="RT Quantiles",main="Stim=A, Response=A")
   abline(0,1)
   plot(q.BB.model,qs.BB,xlab="Diffusion Quantiles",ylab="RT Quantiles",main="Stim=B, Response=B")
   abline(0,1)
   plot(q.AB.model,qs.AB,xlab="Diffusion Quantiles",ylab="RT Quantiles",main="Stim=A, Response=B")
   abline(0,1)
   plot(q.BA.model,qs.BA,xlab="Diffusion Quantiles",ylab="RT Quantiles",main="Stim=B, Response=A")
   abline(0,1)   