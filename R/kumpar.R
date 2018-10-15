# Kumarasuammy Pareto ---------------------------
# Density KumPar
dkumpar<-function(x, beta, k, lambda, phi, log = FALSE)
{
  dens <- (lambda*phi*k*beta^k)/(x^(k+1))*(1-(beta/x)^k)^(lambda-1)*(1-(1-(beta/x)^k)^lambda)^(phi-1)
  if(log[1]) dens<-log(dens)
  ret <- dens
  ret
}

# Probability KumPar
pkumpar<-function(q, beta, k, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  x<-q
  prob <- 1-(1-(1-(beta/x)^k)^lambda)^phi
  if(!lower.tail[1]) prob <- 1-prob
  if(log.p[1]) prob <- log(prob)
  ret <- prob
  ret
}

# Quantile KumPar
qkumpar<-function(p, beta, k, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  q <- beta/((1-(1-(1-p)^(1/phi))^(1/lambda))^(1/k))
  ret <- q
  ret
}

# Hazard KumPar
hkumpar<-function(q, beta, k, lambda, phi)
{
  d <- (lambda*phi*k*beta^k)/(q^(k+1))*(1-(beta/q)^k)^(lambda-1)*(1-(1-(beta/q)^k)^lambda)^(phi-1)
  p <- 1-(1-(1-(beta/q)^k)^lambda)^phi
  h<-d/p
  ret <- h
  ret
}

# Cumulative Hazard KumPar
#Hkumgg<-function(q, tau, alpha=1, k, lambda, phi)
#{
#  p<-pkumgg(q = q, tau = tau, alpha = alpha, k = k, lambda = lambda, phi = phi, lower.tail=FALSE)
#  -log(p)
#}

# Random KumPar
rkumpar<-function(n, beta, k, lambda, phi, cens.prop=0)
{
  if(cens.prop<0 | cens.prop>1){stop("cens.prop must be a numeric between 0 and 1.")}

  u <- runif(n)
  q <- beta/((1-(1-(1-u)^(1/phi))^(1/lambda))^(1/k))
  ret <- q
  if(cens.prop>0){
    quantile<-ret
    ic<-sample(1:n,trunc(n*cens.prop))
    cens.ind<-rep(0,n)
    cens.ind[ic]<-1
    ret<-as.matrix(cbind(quantile,cens.ind))
  }
  ret
}

# Maximum Likelihood Parameter estimation KumPar
