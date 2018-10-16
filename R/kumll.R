#==========================================================
# Kumarasuammy Log-Logistic
#==========================================================

#' Kumarasuammy Log-Logistic Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the KumLL distribution with parameters \code{lambda}, \code{phi},
#' \code{c}, \code{k} and \code{s}.
#'
#' The KumLL distribution was described by Santana \emph{et al} (2012) and has density
#'
#' \deqn{f(x) = (\lambda\phi\gamma)/(\alpha^(\lambda\gamma))x^(\lambda\gamma-1)
#' (1+(t/\alpha)^\gamma)^(-\lambda-1)(1-(1-1/(1+(t/\alpha)^\gamma))^\lambda)^(\phi-1)}
#' with scale parameter \eqn{\alpha}, shape parameters \eqn{\lambda}, \eqn{\phi} and
#' \eqn{\gamma} that govern the distribution's skewness. The parameters \eqn{\lambda}
#' and \eqn{phi}, come from the Kumaraswamy Generalized family introduced by Cordeiro
#' and Castro (2011).
#'
#' The KumLL is a special case of KumBII introduced by Parnaíba \emph{et al} (2013).
#'
#' With \code{phi = 1} KumLL becomes the Exponentiated Log-Logistic distribution.
#' In addition, when \code{lambda = 1} it becomes the Log-Logistic distribution.
#' Those are arguably the most important sub-models to KumLL.
#'
#' When \code{lambda = 1} then the KumLL distribution becomes the BXII distribution
#' described by Zimmer \emph{et al} (1998).
#'
#' This distribution's failure rate function accommodates increasing, decreasing, unimodal
#' and bathtub shaped forms, that depend basically on the values of the shape parameters.
#' Moreover, it is quite flexible for modeling survival data.
#'
#' @return \code{dkumll} gives the density, \code{pkumll} gives the distribution
#' function, \code{qkumll} gives the quantile function, and \code{rkumll} generates random values.
#'
#' The length of the result is determined by \code{n} for \code{rkumll}, for the other fucntions the
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
#' DE SANTANA, T. V. F.; Ortega, E. M.; Cordeiro, G. M.;  Silva, G. O. The
#' Kumaraswamy-log-logistic distribution. Journal of Statistical Theory
#' and Applications, 2012, 11.3: 265-291.
#'
#' PARANAÍBA, P. F.; Ortega, E. M.; Cordeiro, G. M.; Pascoa, M. A. D. The Kumaraswamy
#' Burr XII distribution: theory and practice. Journal of Statistical Computation
#' and Simulation, 2013, 83.11: 2117-2143.
#'
#' CORDEIRO, G. M.; DE CASTRO, M. A new family of generalized distributions. Journal of
#' statistical computation and simulation, 2011, 81.7: 883-898.
#'
#' ZIMMER, W. J.; KEATS, J. B.; WANG, F. K. The Burr XII distribution in reliability
#' analysis. Journal of quality technology, 1998, 30.4: 386-394.
#'
#' @seealso LINK TO OTHER PACKAGE DISTRIBUTIONS
#'
#' @examples
#' # Generating values and comparing with the function
#' x <- rkumll(10000, alpha = 0.5, gamma = 2, lambda = 2, phi = 2)
#' hist(x, probability = T, breaks = 100)
#' curve(dkumll(x, alpha = 0.5, gamma = 2, lambda = 2, phi = 2),
#'       from = 0, to = 25, add = T)
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
#' @name KumLL
NULL

# Density KumLL
#' @rdname KumLL
dkumll<-function(x, alpha, gamma, lambda, phi, log = FALSE)
{
  dens <- (lambda*phi*gamma)/alpha^(lambda*gamma)*x^(lambda*gamma-1)*(1+(x/alpha)^gamma)^(-(lambda+1))*(1-(1-1/(1+(x/alpha)^gamma))^lambda)^(phi-1)
  if(log[1]) dens<-log(dens)
  ret <- dens
  ret
}

# Probability KumLL
#' @rdname KumLL
pkumll<-function(q, alpha, gamma, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  x<-q
  prob <- 1-(1-(1-1/(1+(x/alpha)^gamma))^lambda)^phi
  if(!lower.tail[1]) prob <- 1-prob
  if(log.p[1]) prob <- log(prob)
  ret <- prob
  ret
}

# Quantile KumLL
#' @rdname KumLL
qkumll<-function(p, alpha, gamma, lambda, phi, lower.tail=TRUE, log.p = FALSE)
{
  q <- alpha*(1/(1-(1-(1-p)^(1/phi))^(1/lambda))-1)^(1/gamma)
  ret <- q
  ret
}

# Hazard KumLL
#' @rdname KumLL
hkumll<-function(q, alpha, gamma, lambda, phi)
{
  d <- (lambda*phi*gamma)/alpha^(lambda*gamma)*q^(lambda*gamma-1)*(1+(q/alpha)^gamma)^(-(lambda+1))*(1-(1-1/(1+(q/alpha)^gamma))^lambda)^(phi-1)
  p <- 1-(1-(1-1/(1+(q/alpha)^gamma))^lambda)^phi
  #h<-(lambda*phi*gamma)/alpha^(lambda*gamma)*q^(lambda*gamma-1)*(1+(q/alpha))^(-(lambda+1)))*(1-(1-(1/(a+(q/alpha)^gamma)))^lambda)^(-1)
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

# Random KumLL
#' @rdname KumLL
rkumll<-function(n, alpha, gamma, lambda, phi, cens.prop=0)
{
  if(cens.prop<0 | cens.prop>1){stop("cens.prop must be a numeric between 0 and 1.")}

  u <- runif(n)
  q <- alpha*(1/(1-(1-(1-u)^(1/phi))^(1/lambda))-1)^(1/gamma)
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

# Maximum Likelihood Parameter estimation KumLL
#' @rdname KumLL
mlkumll <- function(x, a.ini, g.ini, l.ini, p.ini){

  # Likelihood Function KumLL
  lkumll<-function(param, x){
    if(any(c(param[1]<0,param[2]<0,param[3]<0,param[4]<0))){
      vero<-NA
    }
    else{
      a<-param[1]
      g<-param[2]
      l<-param[3]
      p<-param[4]

      if(any(param < 1e-20)) return(.Machine$double.xmax^.5)

      f <- dkumll(x, a, g, l, p)

      loglike <- sum(-log(f))
    }
    return(loglike)
  }

  # Estimating parameters
  ini.val <- c(a = 1.5, g = 5.5, l = 1.8, p = 1.1)
  estim <- optim(fn = lkumll, par = ini.val, x = x, hessian = T)
  return(estim)
}
