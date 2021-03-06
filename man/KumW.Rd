% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kumw.R
\name{KumW}
\alias{KumW}
\alias{dkumw}
\alias{pkumw}
\alias{qkumw}
\alias{hkumw}
\alias{rkumw}
\alias{mlkumw}
\title{Kumarasuammy Weibull Distribution}
\source{
The source code of all distributions in this package can also be
found on the \href{https://github.com/aneisse/survdistr}{survdistr} Github repository.
}
\usage{
dkumw(x, beta, c, lambda, phi, log = FALSE)

pkumw(q, beta, c, lambda, phi, lower.tail = TRUE, log.p = FALSE)

qkumw(p, beta, c, lambda, phi, lower.tail = TRUE, log.p = FALSE)

hkumw(q, beta, c, lambda, phi)

rkumw(n, beta, c, lambda, phi, cens.prop = 0)

mlkumw(x, b.ini, c.ini, l.ini, p.ini)
}
\arguments{
\item{x, q}{numeric vector of quantiles.}

\item{beta}{scale parameter \eqn{\beta > 0}.}

\item{c}{shape parameter \eqn{c > 0}.}

\item{lambda}{shape parameter \eqn{\lambda > 0}.}

\item{phi}{shape parameter \eqn{\phi \ge 0}.}

\item{log, log.p}{logical; if \code{TRUE}, probabilities/densities \code{p} are given as \code{log(p)}.}

\item{lower.tail}{logical; if \code{TRUE}, probabilities are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}}

\item{n}{desired size of the random number sample.}

\item{cens.prop}{proportion of censored data to be simulated. If greater than \code{0}, a matrix will be
returned instead of a vector. The matrix will contain the random values and a censorship indicator variable.}
}
\value{
\code{dkumw} gives the density, \code{pkumw} gives the distribution
function, \code{qkumw} gives the quantile function, and \code{rkumw} generates random values.

The length of the result is determined by \code{n} for \code{rkumw}, for the other fucntions the
length is the same as the vector passed to the first argument.

Only the first element of the logical arguments are used.
}
\description{
Density, distribution function, quantile function and random generation
for the KumW distribution with parameters \code{lambda}, \code{phi},
\code{c}, \code{k} and \code{s}.
}
\details{
The KumW distribution was described by Cordeiro \emph{et al} (2010) and has density

\deqn{f(x) = \lambda\phic\beta^cx^(x-1)e^(-(\betax)^c)(1-e^(-(\betax)^c))^(\lambda-1)
(1-(1-e^(-(\betax)^c))^\lambda)^(\phi-1)}
with scale parameter \eqn{\beta}, shape parameters \eqn{\lambda}, \eqn{\phi} and
\eqn{c}. The parameters \eqn{\lambda} and \eqn{phi}, come from the Kumaraswamy
Generalized family introduced by Cordeiro and Castro (2011).

The KumW is a special case of KumBII introduced by Parnaíba \emph{et al} (2013).

With \code{phi = 1} KumW becomes the Exponentiated Weibull distribution.
In addition, when \code{lambda = 1} it becomes the Weibull distribution.

When \code{phi = c = 1} then the KumW distribution becomes the exponentiated exponential
distribution .

The above are arguably the most important sub-models to KumW. More su-models are decribed
by Cordeiro \emph{et al} (2010) as well as some expasions for the KumW pdf.
}
\examples{
# Equivalency with the Weibull
all.equal(dkumw(5, beta = 0.5, c = 2, lambda = 1, phi = 1),
          dweibull(5, shape = 2, scale = 1/0.5))

# Generating values and comparing with the function
x <- rkumw(10000, beta = 1.5, c = 0.5, lambda = 3, phi = 10)
hist(x, probability = T, breaks = 100)
curve(dkumw(x, beta = 1.5, c = 0.5, lambda = 3, phi = 10),
      from = 0, to = 2, add = T)

}
\references{
CORDEIRO, G. M.; ORTEGA, E. M. M; NADARAJAH, S.. The Kumaraswamy Weibull
distribution with application to failure data. Journal of the Franklin
Institute, 2010, 347.8: 1399-1429.

PARANAÍBA, P. F.; Ortega, E. M.; Cordeiro, G. M.; Pascoa, M. A. D. The Kumaraswamy
Burr XII distribution: theory and practice. Journal of Statistical Computation
and Simulation, 2013, 83.11: 2117-2143.

CORDEIRO, G. M.; DE CASTRO, M. A new family of generalized distributions. Journal of
statistical computation and simulation, 2011, 81.7: 883-898.
}
\seealso{
LINK TO OTHER PACKAGE DISTRIBUTIONS
}
\author{
Anderson Neisse <a.neisse@gmail.com>
}
