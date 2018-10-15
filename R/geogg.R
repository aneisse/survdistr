#==========================================================
# Geometrical Generalized Gamma
#==========================================================

#' Geometrical Generalized Gamma Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the GeoGG distribution with parameters \code{tau}, \code{alpha},
#' \code{k} and \code{pg}.
#'
#' The GeoGG distribution has density
#'
#' \deqn{f(x) = ((\tau(1-p))/(\alpha \Gamma(k)))(x/\alpha))^(\tau k-1)
#' e^(-(t/\alpha)^\tau)(1 - p*(1-\gamma(k, (t/\alpha)^\tau)))^(-2)}
#' with scale parameter \eqn{\alpha}, shape parameters \eqn{\tau} and \eqn{k} and \eqn{p}
#' is the degeneration parameter as described by Ortega \emph{et al} (2011).
#'
#' If \code{alpha} is not specified it assumes the default value of 1.
#'
#' With \code{pg = 0} GeoGG equals the Generalized Gamma with parameterizantion
#' described on Stacy's (1962).
#'
#' When \code{pg = 0} and \code{tau = 1} then it equals a Gama
#' with \code{shape = k} and \code{scale = alpha}.
#'
#' With \code{pg = 0} and \code{k = 1} it equals a Weibull with
#' \code{shape = tau} and \code{scale = alpha}.
#'
#' Finally, when \code{pg = 0} and \code{tau = k = 1} it becomes the
#' Exponential with \code{rate = 1/alpha}.
#'
#' For the cases described above, when \code{0 < pg < 1} then the GeoGG will
#' result in the Geometrical Gamma, Geometrical Weibull and Geometrical Exponential,
#' respectively, as described by Ortega \emph{et al} (2011).
#'
#' @return \code{dgeogg} gives the density, \code{pgeogg} gives the distribution
#' function, \code{qgeogg} gives the quantile function, and \code{rgeogg} generates random values.
#'
#' The length of the result is determined by \code{n} for \code{rgeogg}, for the other fucntions the
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
#' STACY, E. W., et al. A generalization of the gamma distribution. The Annals
#' of mathematical statistics, 1962, 33.3: 1187-1192.
#'
#' ORTEGA, E. M. M.; CORDEIRO, G. M.; PASCOA, M. A. R. The generalized gamma
#' geometric distribution. Journal of Statistical Theory and Applications,
#' 2011, 10.3: 433-454.
#'
#' @seealso LINK TO OTHER PACKAGE DISTRIBUTIONS
#'
#' @examples
#' #Equivalency to the Generalized Gamma on it's original parameterization
#' all.equal(dgeogg(5, alpha = 2, tau = 3, k = 0.5, pg = 0),
#'           dgengamma.orig(5, shape = 3, scale = 2, k = 0.5))
#'
#' #Equivalency to the Gamma
#' all.equal(dgeogg(5, alpha = 2, tau = 1, k = 0.5, pg = 0),
#'           dgamma(5, shape = 0.5, scale = 2))
#'
#' #Equivalency to the Weibull
#' all.equal(dgeogg(5, alpha = 2, tau = 3, k = 1, pg = 0),
#'           dweibull(5, shape = 3, scale = 2))
#'
#' #Equivalency to the Exponential
#' all.equal(dgeogg(5, alpha = 2, tau = 1, k = 1, pg = 0),
#'           dexp(5, rate = 1/2))
#'
#' # Generating values and comparing with the function
#' x <- rgeogg(10000, alpha = 2, tau = 3, k = 0.5, pg = 0.35)
#' hist(x, probability = T, breaks = 100)
#' curve(dgeogg(x, alpha = 2, tau = 3, k = 0.5, pg = 0.35),
#'       from = 0, to = 4, add = T)
#'
#' @param x,q numeric vector of quantiles.
#' @param alpha scale parameter \eqn{\alpha} from the Generalized Gama (Stacy, 1962), \eqn{\alpha > 0}.
#' @param tau shape parameter \eqn{\tau} from the Generalized Gama (Stacy, 1962), \eqn{\tau > 0}.
#' @param k shape parameter \eqn{k} from the Generalized Gama (Stacy, 1962), \eqn{k > 0}.
#' @param pg is the probability of success \eqn{p} from the Geometric, \eqn{0 < p < 1}.
#' @param n desired size of the random number sample.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \code{p} are given as \code{log(p)}.
#' @param cens.prop proportion of censored data to be simulated. If greater than \code{0}, a matrix will be
#' returned instead of a vector. The matrix will contain the random values and a censorship indicator variable.
#' @name GeoGG
NULL

# Density GeoGG
#' @rdname GeoGG
dgeogg<-function(x, alpha=1, tau, k, pg, log = FALSE)
{
  gr <- pgengamma.orig(q = x, shape = tau, scale = alpha, k = k)
  dens <-((tau*(1-pg))/(alpha*gamma(k))) * (x/alpha)^(tau*k-1) * exp(-(x/alpha)^tau) * (1-pg*(1-gr))^(-2)
  if(log[1]) dens<-log(dens)
  ret <- dens
  ret
}

# Probability GeoGG
#' @rdname GeoGG
pgeogg<-function(q, alpha=1, tau, k, pg, lower.tail=TRUE, log.p = FALSE)
{
  gr <- pgengamma.orig(q = q, shape = tau, scale = alpha, k = k)
  prob <- gr/(1-pg*(1-gr))
  if(!lower.tail[1]) prob <- 1-prob
  if(log.p[1]) prob <- log(prob)
  ret <- prob
  ret
}

# Quantile GeoGG
#' @rdname GeoGG
qgeogg<-function(p, alpha=1, tau, k, pg, lower.tail=TRUE, log.p = FALSE)
{
  u <- (p*pg-p)/(p*pg-1)
  q <- qgengamma.orig(p=u, shape = tau, scale = alpha, k = k, lower.tail = lower.tail, log.p = log.p)
  ret <- q
  ret
}

# Random GeoGG
#' @rdname GeoGG
rgeogg<-function(n, alpha=1, tau, k, pg,cens.prop=0)
{
  if(cens.prop<0 | cens.prop>1){stop("cens.prop must be a numeric between 0 and 1.")}

  p <- runif(n)
  u <- (p*pg-p)/(p*pg-1)
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

# Maximum Likelihood Parameter estimation GeoGG
#' @rdname GeoGG
ml.geogg<-function(x, alpha.ini, tau.ini, k.ini, pg.ini)
{
  l.geogg<-function(param)
  {
    tau <- param[1]
    alpha <- param[2]
    k <- param[3]
    pg <- param[4]
    if(any(tau < 1e-20)) return(.Machine$double.xmax^.5)
    if(any(alpha < 1e-20)) return(.Machine$double.xmax^.5)
    if(any(k < 1e-20)) return(.Machine$double.xmax^.5)
    if(any(pg < 1e-20)) return(.Machine$double.xmax^.5)
    gr <- pgengamma.orig(q = x, shape = tau, scale = alpha, k = k)
    f <- ((tau*(1-pg))/(alpha*gamma(k))) * (x/alpha)^(tau*k-1) * exp(-(x/alpha)^tau) * (1-pg*(1-gr))^(-2)
    f <- f[f > 0 & f < Inf]
    f <- log(f + 2.225074e-308)
    v <- sum(-f, na.rm = TRUE)
    #if(any("parametro errado")){v <- NA}
  }

  ml <- optim(par = c(tau.ini,alpha.ini,k.ini,pg.ini), fn = l.geogg ,gr = NULL, method="BFGS",hessian=T,control = list(type = 1))
  t<-list("par"=ml$par,"Convergence"=ml$convergence,"Hessian"=ml$hessian)
  return(t)
}
