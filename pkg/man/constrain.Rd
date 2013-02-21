\name{constrain}
\alias{constrain}
\title{constrain}
\description{Make a constraint list. Write each constraint as c_1 v_1 + ... +
c_n v_n (in)equality constant, where c_i are real coefficients and
v_i are optimization variables.}
\usage{constrain(equality, constant, ...)}
\arguments{
  \item{equality}{"eq" or "geq" designating whether this is an equality or an
inequality constraint.}
  \item{constant}{constant on the right side of the inequality.}
  \item{\dots}{Variables to include in this constraint, usually using something
like variable("xi",i,j,-1).}
}



\author{Toby Dylan Hocking}




\examples{
constraints <- list()
for(i in 1:3){
  constraints <-
    c(constraints,list(constrain("geq",0,variable("xi",i))))
}

lapply(1:3,function(i)constrain("geq",0,variable("xi",i)))
}
