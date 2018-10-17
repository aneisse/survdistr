#==========================================================
# Kumarasuammy Weibull
#==========================================================

#' Kumarasuammy Weibull Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the KumW distribution with parameters \code{lambda}, \code{phi},
#' \code{c}, \code{k} and \code{s}.
#'
#' The KumW distribution was described by Cordeiro \emph{et al} (2010) and has density
#'
#' \deqn{f(x) = \lambda\phic\beta^cx^(x-1)e^(-(\betax)^c)(1-e^(-(\betax)^c))^(\lambda-1)
#' (1-(1-e^(-(\betax)^c))^\lambda)^(\phi-1)}
#' with scale parameter \eqn{\beta}, shape parameters \eqn{\lambda}, \eqn{\phi} and
#' \eqn{c}. The parameters \eqn{\lambda} and \eqn{phi}, come from the Kumaraswamy
#' Generalized family introduced by Cordeiro and Castro (2011).
#'
#' The KumW is a special case of KumBII introduced by Parnaíba \emph{et al} (2013).
#'
#' With \code{phi = 1} KumW becomes the Exponentiated Weibull distribution.
#' In addition, when \code{lambda = 1} it becomes the Weibull distribution.
#'
#' When \code{phi = c = 1} then the KumW distribution becomes the exponentiated exponential
#' distribution .
#'
#' The above are arguably the most important sub-models to KumW. More su-models are decribed
#' by Cordeiro \emph{et al} (2010) as well as some expasions for the KumW pdf.
#'
#' @return \code{dkumw} gives the density, \code{pkumw} gives the distribution
#' function, \code{qkumw} gives the quantile function, and \code{rkumw} generates random values.
#'
#' The length of the result is determined by \code{n} for \code{rkumw}, for the other fucntions the
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
#' CORDEIRO, G. M.; ORTEGA, E. M. M; NADARAJAH, S.. The Kumaraswamy Weibull
#' distribution with application to failure data. Journal of the Franklin
#' Institute, 2010, 347.8: 1399-1429.
#'
#' PARANAÍBA, P. F.; Ortega, E. M.; Cordeiro, G. M.; Pascoa, M. A. D. The Kumaraswamy
#' Burr XII distribution: theory and practice. Journal of Statistical Computation
#' and Simulation, 2013, 83.11: 2117-2143.
#'
#' CORDEIRO, G. M.; DE CASTRO, M. A new family of generalized distributions. Journal of
#' statistical computation and simulation, 2011, 81.7: 883-898.
#'
#' @seealso LINK TO OTHER PACKAGE DISTRIBUTIONS
#'
#' @examples
#' # Equivalency with the Weibull
#' all.equal(dkumw(5, beta = 0.5, c = 2, lambda = 1, phi = 1),
#'           dweibull(5, shape = 2, scale = 1/0.5))
#'
#' # Generating values and comparing with the function
#' x <- rkumw(10000, beta = 1.5, c = 0.5, lambda = 3, phi = 10)
#' hist(x, probability = T, breaks = 100)
#' curve(dkumw(x, beta = 1.5, c = 0.5, lambda = 3, phi = 10),
#'       from = 0, to = 2, add = T)
#'
#' @param x,q numeric vector of quantiles.
#' @param lambda shape parameter \eqn{\lambda > 0}.
#' @param phi shape parameter \eqn{\phi \ge 0}.
#' @param beta scale parameter \eqn{\beta > 0}.
#' @param c shape parameter \eqn{c > 0}.
#' @param n desired size of the random number sample.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \code{p} are given as \code{log(p)}.
#' @param cens.prop proportion of censored data to be simulated. If greater than \code{0}, a matrix will be
#' returned instead of a vector. The matrix will contain the random values and a censorship indicator variable.
#' @name KumW
NULL

# Density Kum Weibull
#' @rdname KumW
dkumw <- function(x, beta, c, lambda, phi, log = FALSE)
{
  dens <- lambda*phi*dweibull(x = x,shape = c,scale = 1/beta)*pweibull(q = x,shape = c,scale = 1/beta)^(lambda-1)*(1-pweibull(q = x,shape = c,scale = 1/beta)^lambda)^(phi-1)
  if(log[1]) dens<-log(dens)
  ret <- dens
  ret
}

# Probability Kum Weibull
#' @rdname KumW
pkumw<-function(q, beta, c, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  x<-q
  prob <- 1-(1-(pweibull(q = x, shape = c, scale = 1/beta))^lambda)^phi
  if(!lower.tail[1]) prob <- 1-prob
  if(log.p[1]) prob <- log(prob)
  ret <- prob
  ret
}

# Quantile Kum Weibull
#' @rdname KumW
qkumw<-function(p, beta, c, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  p <- (1-(1-p)^(1/phi))^(1/lambda)
  q <- qweibull(p = p, shape = c, scale = 1/beta)
  ret <- q
  ret
}

# Hazard Kum Weibull
#' @rdname KumW
hkumw<-function(q, beta, c, lambda, phi)
{
  d <- lambda*phi*c*beta^c*q^(c-1)*exp(-(beta*q)^c)*(1-exp(-(beta*q)^c))^(lambda-1)*(1-(1-exp(-(beta*q)^c))^lambda)^(phi-1)
  p <- 1-(1-(1-exp(-(beta*q)^c))^lambda)^phi
  h<-d/p
  ret <- h
  ret
}

# Cumulative Hazard Kum Weibull
#' @rdname KumW
#Hkumgg<-function(q, tau, alpha=1, k, lambda, phi)
#{
#  p<-pkumgg(q = q, tau = tau, alpha = alpha, k = k, lambda = lambda, phi = phi, lower.tail=FALSE)
#  -log(p)
#}

# Random Kum Weibull
#' @rdname KumW
rkumw<-function(n, beta, c, lambda, phi, cens.prop=0)
{
  if(cens.prop<0 | cens.prop>1){stop("cens.prop must be a numeric between 0 and 1.")}

  u <- runif(n)
  p <- (1-(1-u)^(1/phi))^(1/lambda)
  q <- qweibull(p = p, shape = c, scale = 1/beta)
  #q<-  (-log(1-(1-(1-u)^(1/phi))^(1/lambda)))^(1/c)/beta
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

# Maximum Likelihood Parameter estimation KumW
#' @rdname KumW
mlkumw <- function(x, b.ini, c.ini, l.ini, p.ini){

  # Likelihood Function KumLL
  lkumw<-function(param, x){
    if(any(param < 0)){
      vero <- NA
    }
    else{
      b <- param[1]
      c <- param[2]
      l <- param[3]
      p <- param[4]

      if(any(param < 1e-20)) return(.Machine$double.xmax^.5)

      f <- dkumw(x, b, c, l, p)

      loglike <- sum(-log(f))
      return(loglike)
    }
  }

  # Estimating parameters
  ini.val <- c(a = b.ini, b = g.ini, l = l.ini, p = p.ini)
  estim <- optim(fn = lkumw, par = ini.val, x = x, hessian = T)
  return(estim)
}

