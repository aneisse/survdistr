% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kumpar.R
\name{KumPar}
\alias{KumPar}
\alias{dkumpar}
\alias{pkumpar}
\alias{qkumpar}
\alias{hkumpar}
\alias{rkumpar}
\title{Kumarasuammy Pareto Distribution}
\source{
The source code of all distributions in this package can also be
found on the \href{https://github.com/aneisse/survdistr}{survdistr} Github repository.
}
\usage{
dkumpar(x, beta, k, lambda, phi, log = FALSE)

pkumpar(q, beta, k, lambda, phi, lower.tail = TRUE, log.p = FALSE)

qkumpar(p, beta, k, lambda, phi, lower.tail = TRUE, log.p = FALSE)

hkumpar(q, beta, k, lambda, phi)

rkumpar(n, beta, k, lambda, phi, cens.prop = 0)
}
\arguments{
\item{x, q}{numeric vector of quantiles \eqn{x > \beta}.}

\item{beta}{scale parameter \eqn{\beta > 0}.}

\item{k}{shape parameter \eqn{\k > 0}.}

\item{lambda}{shape parameter \eqn{\lambda > 0}.}

\item{phi}{shape parameter \eqn{\phi \ge 0}.}

\item{log, log.p}{logical; if \code{TRUE}, probabilities/densities \code{p} are given as \code{log(p)}.}

\item{lower.tail}{logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}}

\item{n}{desired size of the random number sample.}

\item{cens.prop}{proportion of censored data to be simulated. If greater than \code{0}, a matrix will be
returned instead of a vector. The matrix will contain the random values and a censorship indicator variable.}
}
\value{
\code{dkumpar} gives the density, \code{pkumpar} gives the distribution
function, \code{qkumpar} gives the quantile function, and \code{rkumpar} generates random values.

The length of the result is determined by \code{n} for \code{rkumpar}, for the other fucntions the
length is the same as the vector passed to the first argument.

Only the first element of the logical arguments are used.
}
\description{
Density, distribution function, quantile function and random generation
for the KumPar distribution with parameters \code{lambda}, \code{phi},
\code{beta} and \code{k}.
}
\details{
The KumLL distribution was described by Pereira \emph{et al} (2012) and has density

\deqn{f(x) = (\lambda\phik\beta^k)/(x^(k+1))(1-(\beta/x)^k)^(\lambda-1)
(1-(1-(\beta/x)^k)^\lambda)^(\phi-1)}
for \eqn{x > \beta} and with scale parameter \eqn{\beta}, shape parameters \eqn{\lambda}, \eqn{\phi} and
\eqn{k}.

The parameters \eqn{\lambda} and \eqn{phi}, come from the Kumaraswamy
Generalized family introduced by Cordeiro and Castro (2011).

With \code{phi = 1} KumPar becomes the Exponentiated Pareto distribution.
In addition, when \code{lambda = 1} it becomes the Pareto distribution.
}
\examples{
# Generating values and comparing with the function
x <- rkumpar(10000, beta = 2, k = 0.5, lambda = 3, phi = 10)
hist(x, probability = T, breaks = 100)
curve(dkumpar(x, beta = 2, k = 0.5, lambda = 3, phi = 10),
      from = 2, to = 80, add = T)

}
\references{
PEREIRA, M. B.; SILVA, R. B.; ZEA, L. M.; CORDEIRO, G. M. The kumaraswamy Pareto
distribution. arXiv preprint arXiv:1204.1389, 2012.

CORDEIRO, G. M.; DE CASTRO, M. A new family of generalized distributions. Journal of
statistical computation and simulation, 2011, 81.7: 883-898.
}
\seealso{
LINK TO OTHER PACKAGE DISTRIBUTIONS
}
\author{
Anderson Neisse <a.neisse@gmail.com>
}
