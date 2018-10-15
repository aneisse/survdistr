#==========================================================
# Generalized Modified Weibull
#==========================================================

#' Generalized Modified Weibull Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the GMWeibull distribution with parameters \code{alpha}, \code{gamma},
#' \code{lambda} and \code{beta}.
#'
#' The GMWeibull distribution was described by Carrasco \emph{et al} (2008) and has density
#'
#' \deqn{f(x) = (\alpha \beta x^(\gamma-1)(\gamma + \lambda x)e^(\lambda x - \alpha x^\gamma
#' e^(\lambda x)))/(1-e^(-\alpha x^\gamma e^(\lambda x)))^(1-\beta)}
#' with scale parameter \eqn{\alpha}, shape parameters \eqn{\gamma} and \eqn{\beta}. The
#' parameter \eqn{\lambda}, according to the authors, is a kind of
#' accelerating factor working as a parameter of fragility in the survival of the individual as time increases.
#'
#' With \code{lambda = 0} and \code{beta = 1} GMWeibull equals the classical two-parameter Weibull distribution.
#' In addition, if \code{gamma = 1} it equals the Exponential distribution and \code{gamma = 2} it equals the
#' Rayleigh distribution.
#'
#' With \code{gamma = 0} and \code{beta = 1} GMWeibull equals the Extreme-value (log-gamma) distribution.
#'
#' When \code{lambda = 0} GMWeibull reduces to the Exponentiated Weibull distribution described by
#' Mudholkar \emph{et al} (1995). Furthermore, also setting \code{gamma = 1} will reduce GMWeibbul
#' to the exponentiated exponential distribution described by Gupta and Kundu (2001).
#'
#' For \code{beta = 1} GMWeibull will reduce to the Modified Weibull dsitribution introduced by
#' Lai \emph{et al} (2003).
#'
#' @return \code{dgmweibull} gives the density, \code{pgmweibull} gives the distribution
#' function, \code{qgmweibull} gives the quantile function, and \code{rgmweibull} generates random values.
#'
#' The length of the result is determined by \code{n} for \code{rgmweibull}, for the other fucntions the
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
#' CARRASCO, J. M. F.; ORTEGA, E. M. M.; CORDEIRO, G. M. A generalized
#' modified Weibull distribution for lifetime modeling. Computational Statistics
#' & Data Analysis, 2008, 53.2: 450-462.
#'
#' MUDHOLKAR, G. S.; SRIVASTAVA, D. K.; FREIMER, M.. The exponentiated Weibull family:
#' A reanalysis of the bus-motor-failure data. Technometrics, 1995, 37.4: 436-445.
#'
#' LAI, C. D.; XIE, M.; MURTHY, D. N. P. A modified Weibull distribution. IEEE
#' Transactions on reliability, 2003, 52.1: 33-37.
#'
#' @seealso LINK TO OTHER PACKAGE DISTRIBUTIONS
#'
#' @examples
#' # Equivalency with the Exponential distribution
#' all.equal(dgmweibull(5, alpha = 0.5, gamma = 1, lambda = 0, beta = 1),
#'           dexp(5, rate = 0.5))
#'
#' # Generating values and comparing with the function
#' x <- rgmweibull(10000, alpha = 0.5, gamma = 3, lambda = 2, beta = 0.2)
#' hist(x, probability = T, breaks = 100)
#' curve(dgmweibull(x, alpha = 0.5, gamma = 3, lambda = 2, beta = 0.2),
#'       from = 0, to = 25, add = T)
#'
#' @param x,q numeric vector of quantiles.
#' @param alpha scale parameter \eqn{\alpha > 0}.
#' @param gamma shape parameter \eqn{\gamma \ge 0}.
#' @param lambda fragility factor parameter \eqn{\lambda \ge 0}.
#' @param beta shape parameter \eqn{\beta > 0}.
#' @param n desired size of the random number sample.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \code{p} are given as \code{log(p)}.
#' @param cens.prop proportion of censored data to be simulated. If greater than \code{0}, a matrix will be
#' returned instead of a vector. The matrix will contain the random values and a censorship indicator variable.
#' @name GMWeibull
NULL

# Density GM Weibull
#' @rdname GMWeibull
dgmweibull<-function(x, alpha, gamma, lambda, beta, log = FALSE)
{
  dens <- (alpha*beta*x^(gamma-1)*(gamma+lambda*x)*exp(lambda*x-alpha*x^gamma*exp(lambda*x)))/((1-exp(-alpha*x^gamma*exp(lambda*x)))^(1-beta))
  if(log[1]) dens<-log(dens)
  ret <- dens
  ret
}

# Probability GM Weibull
#' @rdname GMWeibull
pgmweibull<-function(q, alpha, gamma, lambda, beta, lower.tail=TRUE, log.p = FALSE)
{
  x<-q
  prob <- (1-exp(-alpha*x^gamma*exp(lambda*x)))^beta
  if(!lower.tail[1]) prob <- 1-prob
  if(log.p[1]) prob <- log(prob)
  ret <- prob
  ret
}

# Quantile GM Weibull
#' @rdname GMWeibull
qgmweibull<-function(p, alpha, gamma, lambda, beta, lower.tail=TRUE, log.p = FALSE)
{
  q<-NULL
  for(i in 1:length(u)){
    r<-function(x) x^gamma*exp(lambda*x)+1/alpha*log(1-p[i]^(1/beta))
    q<-  c(q,uniroot(f = r,interval = c(0,999999))$root)
  }
  ret <- q
  ret
}

# Random GM Weibull
#' @rdname GMWeibull
rgmweibull<-function(n, alpha, gamma, lambda, beta, cens.prop=0)
{
  if(cens.prop<0 | cens.prop>1){stop("cens.prop must be a numeric between 0 and 1.")}

  u <- runif(n)
  q<-NULL
  for(i in 1:length(u)){
    r <- function(x) x^gamma*exp(lambda*x)+1/alpha*log(1 - u[i]^(1/beta))
    q <- c(q,uniroot(f = r,interval = c(0,999999))$root)
  }
  ret <- q
  if(cens.prop>0){
    quantile <- ret
    ic <- sample(1:n,trunc(n*cens.prop))
    cens.ind <- rep(0,n)
    cens.ind[ic] <- 1
    ret <- as.matrix(cbind(quantile,cens.ind))
  }
  ret
}

# Maximum Likelihood Parameter estimation GM Weibull
