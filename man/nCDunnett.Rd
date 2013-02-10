\name{nCDunnett}
\alias{dNCDun}
\alias{pNCDun}
\alias{qNCDun}
\alias{rNCDun}
\docType{package}
\title{Noncentral Dunnett's test distribution}

\description{
Computes the probability density, the cumulative distribution function and the quantile function and generates random samples for the  noncentral Dunnett's test distribution with the numbers means equal to dimension of the vector \code{rho} that correspond to the correlation structure, the degrees of freedom equal to \code{nu} and the number of points of the Gauss-Legendre quadrature equal to \code{n}. The parameter of noncentrality is the vector \code{delta}, with the same dimension of \code{rho}. The two sided test is choose with the \code{two.sided} option, where the default is true. 
}

\usage{
dNCDun(q, nu, rho=(rep(0.5,times=3)), delta=(rep(0, times=3)), n=32, two.sided = TRUE)
pNCDun(q, nu, rho=(rep(0.5,times=3)), delta=(rep(0, times=3)), n=32, two.sided = TRUE)
qNCDun(p, nu, rho=(rep(0.5,times=3)), delta=(rep(0, times=3)), n=32, two.sided = TRUE)
rNCDun(N = 1, nu = Inf, rho=(rep(0.5,times=3)), delta=(rep(0, times=3)), two.sided=TRUE)
}

\arguments{
  \item{q}{
vector of quantiles \eqn{q \in R}.
}
  \item{p}{
vector of probabilities \eqn{(0, 1)}.
}
  \item{rho}{
  vector of correlations, with length equal or greater than \eqn{1}.
}
 \item{N}{
vector size to be simulated, with \eqn{N > 1}.
}
 \item{nu}{
degrees of freedom \eqn{\nu > 0}.
}
 \item{n}{
number of points of the gaussian quadrature \eqn{n > 2}. 
}     
\item{delta}{
vector of noncentrality parameter. Must be of the same length of \code{rho}. 
} 
\item{two.sided}{if true (default) the two-sided distribution is considered, else the one-sided distribution is considered. 
} 
}



\details{
Assumes n = 32 as default value for dNCDun,  pNCDun and qNCDun. The \code{nu} can be a finite real number or an infinity real number. The default value of \code{nu} is \code{Inf} in \code{rNCDun}. When \code{nu=1}, the convergence of the routines requires \code{n>200} points in the Gaussian quadrature to obtain the desired result  accurately. The cumulative distribution  function of the noncentral unilateral Dunnett's test statistic with finite degrees of freedom \eqn{\nu}  is
\deqn{F(q; r, \nu, \bm{\rho}, \bm{\delta})= \displaystyle \int_0^\infty \int_{-\infty}^\infty \phi(y) \prod_{j=1}^r \Phi\left(\frac{\sqrt{\rho_j} y +  x q-\delta_j}{\sqrt{1-\rho_j}}\right)  f(x;\nu)dy dx,}
where \eqn{\bm{\rho}=[\rho_1, \rho_2, \ldots, \rho_r]^\top} is the correlation vector, \eqn{\bm{\delta}=[\delta_1, \delta_2, \ldots, \delta_r]^\top} is the vector of noncentrality parameter, \eqn{q} is the quantile of  the noncentral unilateral Dunnett's test distribution, \eqn{r} is the numbers of means (or sample size) and \eqn{\nu} is the degrees of freedom of a independent chi-square variable in the studentized process.The \eqn{f(x;\nu)} probability density function is given by
\deqn{f(x; \nu)= \frac{ \nu^{\nu/2} }{\Gamma(\nu/2)2^{\nu/2-1}} x^{\nu-1} e^{-\nu x^2/2}, \quad x \ge 0.}

The cumulative distribution  function of the noncentral unilateral Dunnett's test statistic with infinity degrees of freedom is
\deqn{F(q; r, \nu=\infty, \bm{\rho}, \bm{\delta})=  \int_{-\infty}^\infty \phi(y) \prod_{j=1}^r \Phi\left(\frac{\sqrt{\rho_j} y + q-\delta_j }{\sqrt{1-\rho_j}}\right) dy.}

The cumulative distribution  function of the noncentral bilateral  Dunnett's test statistic with finite degrees of freedom \eqn{\nu} is
\deqn{F(q; r, \nu, \bm{\rho}, \bm{\delta})= \int_0^\infty \int_{-\infty}^\infty \phi(y) \prod_{j=1}^r \left[\Phi\left(\frac{\sqrt{\rho_j} y +  x q-\delta_j}{\sqrt{1-\rho_j}}\right) - \Phi\left(\frac{\sqrt{\rho_j} y -  x q-\delta_j}{\sqrt{1-\rho_j}}\right)\right]  f(x;\nu)dy dx.}

Finally, the cumulative distribution  function of the noncentral bilateral  Dunnett's test statistic with infinity degrees of freedom is
\deqn{F(q; r, \nu=\infty, \bm{\rho}, \bm{\delta})= \int_{-\infty}^\infty \phi(y) \prod_{j=1}^r \left[\Phi\left(\frac{\sqrt{\rho_j} y + q -\delta_j}{\sqrt{1-\rho_j}}\right) - \Phi\left(\frac{\sqrt{\rho_j} y - q-\delta_j }{\sqrt{1-\rho_j}}\right)\right] dy.}
}

\value{\code{dNCDun} gives the density, \code{pNCDun} gives the cumulative distribution function, \code{qNCDun} gives the quantile function, and \code{rNCDun} generates random deviates.}

\references{
Broch, SC; Ferreira, DF. The  noncentral distribution of the Dunnett's statistics revisited. Submmited for publications. 2013.
}

\author{
Siomara Cristina Broch  and Daniel Furtado Ferreira.

Maintainer: Daniel Furtado Ferreira \code{<danielff@dex.ufla.br>}
}

\seealso{
Package homepage: \code{<www.dex.ufla.br/~danielff/r_resources.html>}
}

\examples{
library(nCDunnett)
n <- 500
nu <- 9
rho <- c(0.5,0.5,0.5)
delta <- c(0,0,0)
x <- rNCDun(n, nu, rho, delta)
hist(x)
q <- 2.30
p <- 0.95
prob <- length(x[x<=q])/n
prob
pNCDun(q, nu, rho, delta, 32, TRUE)
dNCDun(q, nu, rho, delta, 32, TRUE)
qNCDun(p, nu, rho, delta, 32, TRUE)

q     <- c(2.34, 4.50, 3.40, 1.77)
p     <- c(0.85, 0.95, 0.34, 0.45)
nu    <- c(Inf, 19, 15, 23)
# unbalanced example
rho   <- c(0.23, 0.25, 0.27, 0.65, 0.34)
delta <- c(0, 0, 0, 0, 0) # central case
pNCDun(q, nu, rho, delta, 32, TRUE)
dNCDun(q, nu, rho, delta, 32, TRUE)
qNCDun(p, nu, rho, delta, 32, TRUE)
}

\keyword{Dunnett}