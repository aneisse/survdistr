#==========================================================
# Kumaraswamy BURR XII
#==========================================================

#' Kumaraswamy BURR XII Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the KumBXII distribution with parameters \code{lambda}, \code{phi},
#' \code{c}, \code{k} and \code{s}.
#'
#' The KumBXII distribution was described by Parnaíba \emph{et al} (2013) and has density
#'
#' \deqn{f(x) = \lambda\phicks^(-c)x^(c-1)(1+(x/s)^c)^(-k-1)(1-(1+(x/s)^c)^(-k))^(\lambda-1)*
#' (1-(1-(1+(x/s)^c)^(-k))^\lambda)^(\phi-1)}
#' with scale parameter \eqn{s}, shape parameters \eqn{\lambda}, \eqn{\phi},
#' \eqn{k} and \eqn{c}. The parameters \eqn{\lambda} and \eqn{phi}, come from
#' the Kumaraswamy Generalized family introduced by Cordeiro and Castro (2011).
#'
#' With \code{lambda = phi = 1} KumBXII becomes the BXII distribution introduced
#' by Zimmer \emph{et al} (1998). For \code{phi = 1} KumBXII equals the
#' Exponentiated BXII distribution.
#'
#' When \code{s = 1/m} and \code{k = 1} KumBXII becomes the Kumaraswamy
#' Log-Logistic (KumLL) dsitribution. Additionally, with \code{lambda = phi = 1}
#' it reduces to the Log-Logistic distribution to the Exponentiated Weibull distribution.
#'
#' For \code{lambda = c = 1} and \code{lambda = phi = c = 1}, it reduces to the Kumaraswamy
#' Pareto type II and Pareto type II distributions, respectively.
#'
#' If \code{k} tends to infinite, it is identical to the Kumaraswamy Weibull (KwW)
#' distribution. In addition, if \code{lambda = phi = 1}, it gives the Weibull distribution.
#'
#' The KwBXII distribution is not only convenient for modelling comfortable unimodal-shaped
#' failure rates, but it is also suitable for testing goodness-of-fit of some special models
#' such as the KwLL, KwW and Weibull distributions
#'
#' @return \code{dkumBXII} gives the density, \code{pkumBXII} gives the distribution
#' function, \code{qkumBXII} gives the quantile function, and \code{rkumBXII} generates random values.
#'
#' The length of the result is determined by \code{n} for \code{rkumBXII}, for the other fucntions the
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
#' PARANAÍBA, P. F.; Ortega, E. M.; Cordeiro, G. M.; Pascoa, M. A. D. The Kumaraswamy
#' Burr XII distribution: theory and practice. Journal of Statistical Computation
#' and Simulation, 2013, 83.11: 2117-2143.
#'
#' CORDEIRO, G. M.; DE CASTRO, M. A new family of generalized distributions. Journal of
#' statistical computation and simulation, 2011, 81.7: 883-898.
#'
#' ZIMMER, W. J.; KEATS, J. B.; WANG, F. K. The Burr XII distribution in
#' reliability analysis. Journal of quality technology, 1998, 30.4: 386-394.
#'
#' @seealso LINK TO OTHER PACKAGE DISTRIBUTIONS
#'
#' @examples
#' # Generating values and comparing with the function
#' x <- rkumBXII(10000, s = 0.5, c = 10, k = 1.5, lambda = 2, phi = 0.2)
#' hist(x, probability = T, breaks = 100)
#' curve(dkumBXII(x, s = 0.5, c = 10, k = 1.5, lambda = 2, phi = 0.2),
#'      from = 0, to = 6, add = T)
#'
#' @param x,q numeric vector of quantiles.
#' @param lambda shape parameter \eqn{\lambda > 0}.
#' @param phi shape parameter \eqn{\phi \ge 0}.
#' @param c shape parameter \eqn{c > 0}.
#' @param k shape parameter \eqn{k > 0}.
#' @param s scale parameter \eqn{s > 0}.
#' @param n desired size of the random number sample.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \code{p} are given as \code{log(p)}.
#' @param cens.prop proportion of censored data to be simulated. If greater than \code{0}, a matrix will be
#' returned instead of a vector. The matrix will contain the random values and a censorship indicator variable.
#' @name KumBXII
NULL

# Density KumBXII
#' @rdname KumBXII
dkumBXII<-function(x, s, c, k, lambda, phi, log = FALSE)
{
  dens <- (lambda*phi*c*k*s^(-c)*x^(c-1))*(1+(x/s)^c)^(-k-1)*(1-(1+(x/s)^c)^(-k))^(lambda-1) *(1-(1-(1+(x/s)^c)^(-k))^lambda)^(phi-1)
  if(log[1]) dens<-log(dens)
  ret <- dens
  ret
}

# Probability KumBXII
#' @rdname KumBXII
pkumBXII<-function(q, s, c, k, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  x<-q
  prob <- 1-(1-(1-(1+(x/s)^c)^(-k))^lambda)^phi
  if(!lower.tail) prob <- 1-prob
  if(log.p[1]) prob <- log(prob)
  ret <- prob
  ret
}

# Quantile KumBXII#
#' @rdname KumBXII
qkumBXII<-function(p, s, c, k, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  q <- s*((1-(1-(1-p)^(1/phi))^(1/lambda))^(-1/k)-1)^(1/c)
  ret <- q
  ret
}

# Random KumBXII
#' @rdname KumBXII
rkumBXII<-function(n, s, c, k, lambda, phi, cens.prop=0)
{
  if(cens.prop<0 | cens.prop>1){stop("cens.prop must be a numeric between 0 and 1.")}

  u <- runif(n)
  q <- s*((1-(1-(1-u)^(1/phi))^(1/lambda))^(-1/k)-1)^(1/c)
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

# Maximum Likelihood Parameter estimation KumBXII
