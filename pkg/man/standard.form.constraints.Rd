\name{standard.form.constraints}
\alias{standard.form.constraints}
\title{standard form constraints}
\description{Convert constraints to standard form matrix and vector.}
\usage{standard.form.constraints(constraints, ids)}
\arguments{
  \item{constraints}{List of constraints created using constrain().}
  \item{ids}{List of optimization variable ids created using make.ids().}
}

\value{List of variable coefficients and constraint constants.}

\author{Toby Dylan Hocking}




\examples{
## example: optimization variable \alpha\in\RR^3 subject to the
## constraint that it must be on the standard simplex:
opt.vars <- make.ids(alpha=3)
constraints <- with(opt.vars,c(alpha >= 0,list(sum(alpha) == 1)))
standard.form.constraints(constraints,opt.vars)

## linear svm example
set.seed(1)
p <- 2
y <- rep(c(-1,1),each=20)
x <- replicate(p,rnorm(length(y),y))
plot(x,col=y+2,asp=1)
n <- nrow(x)

vars <- make.ids(slack=n,intercept=1,normal=p)
constraints <- vars$slack >= 0
for(i in 1:n){
  ivars <- with(vars,intercept*y[i] + sum(normal)*(x[i,]*y[i]) + slack[i])
  constraints <- c(constraints,list(ivars >= 1))
}
solver.args <- standard.form.constraints(constraints,vars)
n.vars <- length(unlist(vars))
Dvec <- rep(1e-6,n.vars)
Dvec[vars$normal] <- 1
D <- diag(Dvec)
d <- rep(0,n.vars)
d[vars$slack] <- -1 ## C == 1

sol <- quadprog::solve.QP(D,d,solver.args$A,solver.args$b0)
slack <- sol$solution[vars$slack]
normal <- sol$solution[vars$normal]
intercept <- sol$solution[vars$intercept]

title(paste("A linear Support Vector Machine (SVM):",
            "margin SVs circled, slack drawn in red for other SVs"))
abline(-intercept/normal[2],-normal[1]/normal[2])
abline((1-intercept)/normal[2],-normal[1]/normal[2],lty="dotted")
abline((-1-intercept)/normal[2],-normal[1]/normal[2],lty="dotted")
f <- function(x)intercept+sum(normal*x)
yfx <- apply(x,1,f)*y
on.margin <- abs(yfx-1)<1e-6
points(x[on.margin,],cex=2,col=y[on.margin]+2)
i <- yfx<1
## these 2 complicated formulas calculate the point on margin where
## the data point starts picking up slack
x1 <- ((y[i]-intercept)*normal[1]-
       normal[1]*normal[2]*x[i,2]+
       normal[2]^2*x[i,1])/(normal[2]^2+normal[1]^2)
x2 <- (y[i]-intercept)/normal[2]-normal[1]/normal[2]*x1
segments(x[i,1],x[i,2],x1,x2,col="red")
l2norm <- function(x)sqrt(sum(x^2))
rbind(apply(cbind(x1,x2)-x[i,],1,l2norm)*l2norm(normal),slack[i])
}
