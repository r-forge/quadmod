\name{make.ids}
\alias{make.ids}
\title{make ids}
\description{Construct matrices and vectors with their unique ids}
\usage{make.ids(...)}
\arguments{
  \item{\dots}{Named lists of integer vectors of length <= 2. Each item
represents an optimization variable and its dimension.}
}

\value{Named list of optimization variable ids, which will be used to
construct the standard form matrices.}

\author{Toby Dylan Hocking}




\examples{
vars <- make.ids(xi=5,beta0=1,beta=2)
with(vars, xi >= 0)
}
