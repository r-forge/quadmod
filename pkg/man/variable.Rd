\name{variable}
\alias{variable}
\title{variable}
\description{Make a list that represents an optimization variable and its
coefficient in the context of a constraint.}
\usage{variable(name, i = 1, j = NULL, coef = 1)}
\arguments{
  \item{name}{Name of the optimization variable, created with make.ids.}
  \item{i}{If the optimization variable is a vector, the element. If it is a
matrix, the row.}
  \item{j}{Column of the optimization variable.}
  \item{coef}{Coefficient of the optimization variable.}
}

\value{List with elements "variable" "i" "j" "coef", of class "variable".}

\author{Toby Dylan Hocking}




\examples{
variable("xi",2,coef=-1)
variable("alpha",3,3)
}
