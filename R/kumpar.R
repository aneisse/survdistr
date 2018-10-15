#==========================================================
# Kumarasuammy Pareto
#==========================================================

#' Kumarasuammy Pareto Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the KumPar distribution with parameters \code{lambda}, \code{phi},
#' \code{beta} and \code{k}.
#'
#' The KumLL distribution was described by Pereira \emph{et al} (2012) and has density
#'
#' \deqn{f(x) = (\lambda\phik\beta^k)/(x^(k+1))(1-(\beta/x)^k)^(\lambda-1)
#' (1-(1-(\beta/x)^k)^\lambda)^(\phi-1)}
#' for \eqn{x > \beta} and with scale parameter \eqn{\beta}, shape parameters \eqn{\lambda}, \eqn{\phi} and
#' \eqn{k}.
#'
#' The parameters \eqn{\lambda} and \eqn{phi}, come from the Kumaraswamy
#' Generalized family introduced by Cordeiro and Castro (2011).
#'
#' With \code{phi = 1} KumPar becomes the Exponentiated Pareto distribution.
#' In addition, when \code{lambda = 1} it becomes the Pareto distribution.
#'
#' @return \code{dkumpar} gives the density, \code{pkumpar} gives the distribution
#' function, \code{qkumpar} gives the quantile function, and \code{rkumpar} generates random values.
#'
#' The length of the result is determined by \code{n} for \code{rkumpar}, for the other fucntions the
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
#' PEREIRA, M. B.; SILVA, R. B.; ZEA, L. M.; CORDEIRO, G. M. The kumaraswamy Pareto
#' distribution. arXiv preprint arXiv:1204.1389, 2012.
#'
#' CORDEIRO, G. M.; DE CASTRO, M. A new family of generalized distributions. Journal of
#' statistical computation and simulation, 2011, 81.7: 883-898.
#'
#' @seealso LINK TO OTHER PACKAGE DISTRIBUTIONS
#'
#' @examples
#' # Generating values and comparing with the function
#' x <- rkumpar(10000, beta = 2, k = 0.5, lambda = 3, phi = 10)
#' hist(x, probability = T, breaks = 100)
#' curve(dkumpar(x, beta = 2, k = 0.5, lambda = 3, phi = 10),
#'       from = 2, to = 80, add = T)
#'
#' @param x,q numeric vector of quantiles \eqn{x > \beta}.
#' @param lambda shape parameter \eqn{\lambda > 0}.
#' @param phi shape parameter \eqn{\phi \ge 0}.
#' @param beta scale parameter \eqn{\beta > 0}.
#' @param k shape parameter \eqn{\k > 0}.
#' @param n desired size of the random number sample.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \code{p} are given as \code{log(p)}.
#' @param cens.prop proportion of censored data to be simulated. If greater than \code{0}, a matrix will be
#' returned instead of a vector. The matrix will contain the random values and a censorship indicator variable.
#' @name KumPar
NULL


# Density KumPar
#' @rdname KumPar
dkumpar<-function(x, beta, k, lambda, phi, log = FALSE)
{
  dens <- (lambda*phi*k*beta^k)/(x^(k+1))*(1-(beta/x)^k)^(lambda-1)*(1-(1-(beta/x)^k)^lambda)^(phi-1)
  if(log[1]) dens<-log(dens)
  ret <- dens
  ret
}

# Probability KumPar
#' @rdname KumPar
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
#' @rdname KumPar
qkumpar<-function(p, beta, k, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  q <- beta/((1-(1-(1-p)^(1/phi))^(1/lambda))^(1/k))
  ret <- q
  ret
}

# Hazard KumPar
#' @rdname KumPar
hkumpar<-function(q, beta, k, lambda, phi)
{
  d <- (lambda*phi*k*beta^k)/(q^(k+1))*(1-(beta/q)^k)^(lambda-1)*(1-(1-(beta/q)^k)^lambda)^(phi-1)
  p <- 1-(1-(1-(beta/q)^k)^lambda)^phi
  h<-d/p
  ret <- h
  ret
}

# Cumulative Hazard KumPar
#' @rdname KumPar
#Hkumgg<-function(q, tau, alpha=1, k, lambda, phi)
#{
#  p<-pkumgg(q = q, tau = tau, alpha = alpha, k = k, lambda = lambda, phi = phi, lower.tail=FALSE)
#  -log(p)
#}

# Random KumPar
#' @rdname KumPar
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
