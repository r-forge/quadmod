\name{standard.form.constraints}
\alias{standard.form.constraints}
\title{standard form constraints}
\description{Convert constraints to standard form matrix and vector.}
\usage{
standard.form.constraints(constraints, ids)
}
\arguments{
  \item{constraints}{List of constraints created using constrain().}
  \item{ids}{List of optimization variable ids created using make.ids().}
}

\value{List of variable coefficients and constraint constants.}

\author{Toby Dylan Hocking <toby.hocking@inria.fr>}




\examples{
## example: optimization variable \alpha\in\RR^3 subject to the
## constraint that it must be on the standard simplex:
opt.vars <- make.ids(alpha=3)
constraints <- list(constrain("geq",0,variable("alpha",1)),
                    constrain("geq",0,variable("alpha",2)),
                    constrain("geq",0,variable("alpha",3)),
                    constrain("eq",1,
                              variable("alpha",1),
                              variable("alpha",2),
                              variable("alpha",3)))
standard.form.constraints(constraints,opt.vars)

## linear svm example
set.seed(1)
p <- 2
y <- rep(c(-1,1),each=20)
x <- replicate(p,rnorm(length(y),y))
plot(x,col=y+2)
n <- nrow(x)

vars <- make.ids(slack=n,intercept=1,normal=p)
constraints <- list()
for(i in 1:n){
  cargs <- list("geq",1,
                variable("slack",i),
                variable("intercept",coef=y[i]))
  for(j in 1:p){
    cargs <- c(cargs,list(variable("normal",j,coef=y[i]*x[i,j])))
  }
  constraints <-
    c(constraints,list(constrain("geq",0,variable("slack",i)),
                       do.call(constrain,cargs)))
}
solver.args <- standard.form.constraints(constraints,vars)
n.vars <- length(unlist(vars))
Dvec <- rep(1e-6,n.vars)
Dvec[vars$normal] <- 1
D <- diag(Dvec)
d <- rep(0,n.vars)
d[vars$slack] <- -1 ## C == -1
library(quadprog)
sol <- solve.QP(D,d,solver.args$A,solver.args$b0)
slack <- sol$solution[vars$slack]
normal <- sol$solution[vars$normal]
intercept <- sol$solution[vars$intercept]

abline(-intercept/normal[2],-normal[1]/normal[2])
abline(1-intercept/normal[2],-normal[1]/normal[2],lty="dotted")
abline(-1-intercept/normal[2],-normal[1]/normal[2],lty="dotted")
f <- function(x)intercept+sum(normal*x)
yfx <- apply(x,1,f)*y
on.margin <- abs(yfx-1)<1e-6
points(x[on.margin,],cex=2,col=y+2)
points(x[yfx<1,],pch=2,col=y+2)
}
