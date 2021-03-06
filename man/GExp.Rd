% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gexp.R
\name{GExp}
\alias{GExp}
\alias{dgexp}
\alias{pgexp}
\alias{qgexp}
\alias{rgexp}
\title{Generalized Exponential Distribution}
\source{
The source code of all distributions in this package can also be
found on the \href{https://github.com/aneisse/survdistr}{survdistr} Github repository.
}
\usage{
dgexp(x, alpha, lambda, mu, log = FALSE)

pgexp(q, alpha, lambda, mu, lower.tail = TRUE, log.p = FALSE)

qgexp(p, alpha, lambda, mu, lower.tail = TRUE, log.p = FALSE)

rgexp(n, alpha, lambda, mu, cens.prop = 0)
}
\arguments{
\item{x, q}{numeric vector of quantiles. \eqn{x > \mu}.}

\item{alpha}{shape parameter \eqn{\alpha > 0}.}

\item{lambda}{shape parameter \eqn{\lambda > 0}.}

\item{mu}{location parameter \eqn{\mu < x}.}

\item{log, log.p}{logical; if \code{TRUE}, probabilities/densities \code{p} are given as \code{log(p)}.}

\item{lower.tail}{logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}}

\item{n}{desired size of the random number sample.}

\item{cens.prop}{proportion of censored data to be simulated. If greater than \code{0}, a matrix will be
returned instead of a vector. The matrix will contain the random values and a censorship indicator variable.}
}
\value{
\code{dgexp} gives the density, \code{pgexp} gives the distribution
function, \code{qgexp} gives the quantile function, and \code{rgexp} generates random values.

The length of the result is determined by \code{n} for \code{rgexp}, for the other fucntions the
length is the same as the vector passed to the first argument.

Only the first element of the logical arguments are used.
}
\description{
Density, distribution function, quantile function and random generation
for the GExp distribution with parameters \code{alpha}, \code{lambda} and \code{mu}.
}
\details{
The GExp distribution has density

\deqn{f(x) = (\alpha/\lambda)(1-e^(-(x-\mu)/\lambda))^(\alpha-1)e^(-(x-\mu)/\lambda)}
with shape parameter \eqn{\alpha}, scale parameter \eqn{\lambda} and location parameter \eqn{\mu}
and \eqn{x > \mu} as described by Gupta and Kundu (1999).

With \code{alpha = 1} GExp equals the Two-parameter Exponential distribution
with \code{rate = 1/lambda}. Such dsitribution can be computed by a Exponential
transforming the variable \code{g(x) = x - mu}.

With \code{alpha = 1} and \code{mu = 0} GExp equals the usual Exponential distribution.
}
\examples{
# Equivalency with the Two-parameter Exponential distribution
all.equal(dgexp(5, alpha = 1, lambda = 2, mu = 3),
          dexp(5 - 3, rate = 1/2))

# Equivalency with the exponential distribution
all.equal(dgexp(5, alpha = 1, lambda = 2, mu = 0),
          dexp(5, rate = 1/2))

# Generating values and comparing with the function
x <- rgexp(10000, alpha = 1.5, lambda = 2, mu = 3)
hist(x, probability = T, ylim = c(0, 0.5), breaks = 100)
curve(dgexp(x, alpha = 1.5, lambda = 2, mu = 3),
      from = 0, to = 25, add = T)

}
\references{
GUPTA, R. D.; KUNDU, D. Theory & methods: Generalized exponential
distributions. Australian & New Zealand Journal of Statistics, 1999,
41.2: 173-188.

LAWLESS, J. F. Prediction intervals for the two parameter exponential
distribution. Technometrics, 1977, 19.4: 469-472.
}
\seealso{
LINK TO OTHER PACKAGE DISTRIBUTIONS
}
\author{
Anderson Neisse <a.neisse@gmail.com>
}
