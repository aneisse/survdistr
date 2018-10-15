#==========================================================
# Kumaraswamy Generalized Gamma
#==========================================================

#' Kumaraswamy Generalized Gamma Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the KumGG distribution with parameters \code{lambda}, \code{phi},
#' \code{tau}, \code{alpha} and \code{k}.
#'
#' The KumGG distribution was described by Ortega \emph{et al} (2011) and has density
#'
#' \deqn{f(x) = (\lambda\phi\tau)/(\alpha\Gamma[k])(x/\alpha)^(\tauk-1)e^(-(x/\alpha)^\tau)
#' (\gamma[k, (x/\alpha)^\tau])^(\lambda-1)(1-(\gamma[k, (x/\alpha)^\tau])^\lambda)^(\phi-1)}
#'
#' where \eqn{\gamma[., .]} is the incomplete gamma ratio and \eqn{\Gamma[.]} is the gamma funcion.
#' The scale parameter is \eqn{\alpha}, the shape parameters are \eqn{\lambda}, \eqn{\phi} and
#' \eqn{\tau} and \eqn{k}. The parameters \eqn{\lambda} and \eqn{phi}, come from the Kumaraswamy
#' Generalized family introduced by Cordeiro and Castro (2011).
#'
#' With \code{phi = 1} KumGG becomes the Exponentiated Generalized Gamma distribution described
#' by Cordeiro \emph{et al} (2011). Additionally, if \code{tau = k = 1} the Exponentiated Exponential.
#'
#' When \code{k = 1} then the KumGG distribution becomes the KumW distribution
#' described by Cordeiro \emph{et al} (2010). For \code{k = 1} the KumGG becomes KumG distribution.
#'
#' The above are arguably the most important sub-models of KumGG. More sub-models are described in
#' Ortega \emph{et al} (2011).
#'
#' @return \code{dkumgg} gives the density, \code{pkumgg} gives the distribution
#' function, \code{qkumgg} gives the quantile function, and \code{rkumgg} generates random values.
#'
#' The length of the result is determined by \code{n} for \code{rkumgg}, for the other fucntions the
#' length is the same as the vector passed to the first argument.
#'
#' Only the first element of the logical arguments are used.
#'
#' @source The source code of all distributions in this package can also be
#' found on the \href{https://github.com/aneisse/survdistr}{survdistr} Github repository.
#'
#' @author Anderson Neisse <a.neisse@@gmail.com>
#'
#' @references
#' ORTEGA, E. M. M.; CORDEIRO, G. M.; PASCOA, M. A. R. The generalized gamma
#' geometric distribution. Journal of Statistical Theory and Applications,
#' 2011, 10.3: 433-454.
#'
#' CORDEIRO, G. M.; ORTEGA, E. M. M; SILVA, G. O. The exponentiated generalized
#' gamma distribution with application to lifetime data. Journal of statistical
#' computation and simulation, 2011, 81.7: 827-842.
#'
#' CORDEIRO, G. M.; ORTEGA, E. M. M; NADARAJAH, S.. The Kumaraswamy Weibull
#' distribution with application to failure data. Journal of the Franklin
#' Institute, 2010, 347.8: 1399-1429.
#'
#' CORDEIRO, G. M.; DE CASTRO, M. A new family of generalized distributions. Journal of
#' statistical computation and simulation, 2011, 81.7: 883-898.
#'
#' @seealso LINK TO OTHER PACKAGE DISTRIBUTIONS
#'
#' @examples
#' # Generating values and comparing with the function
#' x <- rkumgg(10000, tau = 0.9, alpha = 2, k = 0.5, lambda = 3, phi = 10)
#' hist(x, probability = T, breaks = 100)
#' curve(dkumgg(x, tau = 0.9, alpha = 2, k = 0.5, lambda = 3, phi = 10),
#'       from = 0, to = 3, add = T)
#'
#' @param x,q numeric vector of quantiles.
#' @param lambda shape parameter \eqn{\lambda > 0}.
#' @param phi shape parameter \eqn{\phi \ge 0}.
#' @param alpha scale parameter \eqn{\alpha > 0}.
#' @param gamma shape parameter \eqn{\gamma > 0}.
#' @param n desired size of the random number sample.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \code{p} are given as \code{log(p)}.
#' @param cens.prop proportion of censored data to be simulated. If greater than \code{0}, a matrix will be
#' returned instead of a vector. The matrix will contain the random values and a censorship indicator variable.
#' @name KumGG
NULL

# Density KumGG
#' @rdname KumGG
dkumgg<-function(x, tau, alpha=1, k, lambda, phi, log = FALSE)
{
  gr <- pgengamma.orig(q = x, shape = tau, scale = alpha, k = k)
  dens <- ((lambda*phi*tau)/(alpha*gamma(k)))*(x/alpha)^(tau*k-1)*exp(-(x/alpha)^tau)*(gr)^(lambda-1)*(1-(gr)^lambda)^(phi-1)
  if(log[1]) dens<-log(dens)
  ret <- dens
  ret
}

# Probability KumGG
#' @rdname KumGG
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
#' @rdname KumGG
qkumgg<-function(p, tau, alpha=1, k, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  u <- (1-(1-p)^(1/phi))^(1/lambda)
  q <- qgengamma.orig(p=u, shape = tau, scale = alpha, k = k, lower.tail = lower.tail, log.p = log.p)
  ret <- q
  ret
}

# Hazard KumGG
#' @rdname KumGG
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
#' @rdname KumGG
Hkumgg<-function(q, tau, alpha=1, k, lambda, phi)
{
  p<-pkumgg(q = q, tau = tau, alpha = alpha, k = k, lambda = lambda, phi = phi, lower.tail=FALSE)
  -log(p)
}

# Random KumGG
#' @rdname KumGG
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
#' @rdname KumGG
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
