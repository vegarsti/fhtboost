% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_to_FHT_loss_function.R
\name{data_to_FHT_loss_function}
\alias{data_to_FHT_loss_function}
\title{Data to loss function function}
\usage{
data_to_FHT_loss_function(X, Z, times, delta)
}
\arguments{
\item{X}{Design matrix for \code{ln y0 = beta^T * X}}

\item{Z}{Design matrix for \code{   mu = gamma^T * Z}}

\item{times}{Survival times (possibly censored)}

\item{delta}{Vector of same length as \code{t}, indicating if observations are actually observed (1) or censored (0)}
}
\value{
\code{optimizable_function} Function taking parameters as a list (containing \code{beta}, \code{gamma}).
}
\description{
Takes design matrices \code{X}, \code{Z} and data \code{t}, {delta}, and returns a function to be optimized.
}
\examples{
R code here showing how your function works
}
\keyword{keywords}
