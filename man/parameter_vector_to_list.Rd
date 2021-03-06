% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_vector_to_list.R
\name{parameter_vector_to_list}
\alias{parameter_vector_to_list}
\title{Parameter vector to list}
\usage{
parameter_vector_to_list(vec, p, d)
}
\arguments{
\item{vec}{The vector in question}

\item{p}{Second dimension of design matrix \code{X}, i.e., size of \code{beta}.}

\item{d}{Second dimension of design matrix \code{Z}, i.e., size of \code{gamma}.}
}
\value{
\code{beta} Parameter vector

\code{gamma} Parameter vector
}
\description{
Takes a parameter vector (from \code{nlm}) and splits it into correct parameter vectors.
}
\examples{
R code here showing how your function works
}
\keyword{keywords}
