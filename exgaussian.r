```{r}
# A function to simulate the ex-Gaussian
# Mean = mu + tau
rexgauss <- function(N,mu,sigma,tau) {
  if(sigma < 1e-10) 
  {f <- mu + tau*rexp(N,1)} else if (tau < 1e-10) { 
    f <- rnorm(N,mu,sigma) } else {
      f<- rnorm(N,mu,sigma) + tau*rexp(N,1) } 
  return(f)
}

# The ex-Gaussian pdf
dexgauss <- function(t,mu,sigma,tau,log=TRUE) {
  if (sigma < 1e-10) {
    f <- dexp(t - mu,1/tau,log=TRUE)} else if (tau < 1e-10) { 
      f <- dnorm(t,mu,sigma) } else {
        f<- -log(tau) + (mu-t)/tau + sigma^2/2/tau^2 +
          log(pnorm((t-mu)/sigma - sigma/tau))}
  if (log==FALSE) f<-exp(f)
  return(f)
}
```
