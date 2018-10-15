#==========================================================
# Kumaraswamy Generalized Gamma
#==========================================================


# Density KumGG
dkumgg<-function(x, tau, alpha=1, k, lambda, phi, log = FALSE)
{
  gr <- pgengamma.orig(q = x, shape = tau, scale = alpha, k = k)
  dens <- ((lambda*phi*tau)/(alpha*gamma(k)))*(x/alpha)^(tau*k-1)*exp(-(x/alpha)^tau)*(gr)^(lambda-1)*(1-(gr)^lambda)^(phi-1)
  if(log[1]) dens<-log(dens)
  ret <- dens
  ret
}

# Probability KumGG
pkumgg<-function(q, tau, alpha=1, k, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  gr <- pgengamma.orig(q = q, shape = tau, scale = alpha, k = k)
  prob <- 1-((1-gr^lambda)^phi)
  if(!lower.tail) prob <- 1-prob
  if(log.p[1]) prob <- log(prob)
  ret <- prob
  ret
}

# Quantile KumGG
qkumgg<-function(p, tau, alpha=1, k, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  u <- (1-(1-p)^(1/phi))^(1/lambda)
  q <- qgengamma.orig(p=u, shape = tau, scale = alpha, k = k, lower.tail = lower.tail, log.p = log.p)
  ret <- q
  ret
}

# Hazard KumGG
hkumgg<-function(q, tau, alpha=1, k, lambda, phi)
{
  #d<-dkumgg(x = q, tau = tau, alpha = alpha, k = k, lambda = lambda, phi = phi)
  #p<-pkumgg(q = q, tau = tau, alpha = alpha, k = k, lambda = lambda, phi = phi, lower.tail=FALSE)
  #d/p
  gr <- pgengamma.orig(q = q, shape = tau, scale = alpha, k = k)
  h <- ((lambda*phi*tau)/(alpha*gamma(k)))*(q/alpha)^(tau*k-1)*exp(-(q/alpha)^tau)*(gr)^(lambda-1)*(1-(gr)^lambda)^(-1)
  h
}

# Cumulative Hazard KumGG
Hkumgg<-function(q, tau, alpha=1, k, lambda, phi)
{
  p<-pkumgg(q = q, tau = tau, alpha = alpha, k = k, lambda = lambda, phi = phi, lower.tail=FALSE)
  -log(p)
}

# Random KumGG
rkumgg<-function(n, tau, alpha=1, k, lambda, phi,cens.prop=0)
{
  if(cens.prop<0 | cens.prop>1){stop("cens.prop must be a numeric between 0 and 1.")}

  p <- runif(n)
  u <- (1-(1-p)^(1/phi))^(1/lambda)
  q <- qgengamma.orig(p=u, shape = tau, scale = alpha, k = k)
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

# Maximum Likelihood Parameter estimation KumGG
ml.kumgg<-function(x, tau.ini, alpha.ini, k.ini, lambda.ini, phi.ini)
{
  l.kumgg<-function(param)
  {
    tau <- param[1]
    alpha <- param[2]
    k <- param[3]
    lambda <- param[4]
    phi<-param[5]
    if(any(tau < 1e-20)) return(.Machine$double.xmax^.5)
    if(any(alpha < 1e-20)) return(.Machine$double.xmax^.5)
    if(any(k < 1e-20)) return(.Machine$double.xmax^.5)
    if(any(lambda < 1e-20)) return(.Machine$double.xmax^.5)
    if(any(phi < 1e-20)) return(.Machine$double.xmax^.5)
    gr <- pgengamma.orig(q = x, shape = tau, scale = alpha, k = k)
    f <- ((lambda*phi*tau)/(alpha*gamma(k))) * ((x/alpha)^(tau*k-1)) * exp(-(x/alpha)^tau) * gr^(lambda-1) * ((1-gr^lambda)^(phi-1))
    f <- f[f > 0 & f < Inf]
    f <- log(f + 2.225074e-308)
    v <- sum(-f, na.rm = TRUE)
    #if(any("parametro errado")){v <- NA}
  }

  ml <- optim(par = c(tau.ini,alpha.ini,k.ini,lambda.ini,phi.ini), fn = l.kumgg ,gr = NULL, method="BFGS",hessian=T,control = list(type = 1))
  t<-list("par"=ml$par,"Convergence"=ml$convergence,"Hessian"=ml$hessian)
  return(t)
}
