% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loss_function_derivative_mu.R
\name{loss_function_derivative_mu}
\alias{loss_function_derivative_mu}
\title{Loss function for \code{mu}}
\usage{
loss_function_derivative_mu(y0, mu, sigma2, times, delta)
}
\arguments{
\item{y0}{Initial level of underlying Wiener process \code{y0}}

\item{mu}{Drift of underlying Wiener process \code{mu}}
}
\value{
gradient The gradient
}
\description{
Derivative of minus log likelihood with respect to
one parameter: \code{mu}
}
\examples{
R code here showing how your function works
}
\keyword{keywords}
