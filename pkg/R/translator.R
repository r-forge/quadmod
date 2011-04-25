make.ids <- structure(function
### Construct matrices and vectors with their unique ids
(...
### Named lists of integer vectors of length <= 2. Each item
### represents an optimization variable and its dimension.
 ){
  vars <- list(...)
  ids <- list()
  i <- 1
  for(N in names(vars)){
    d <- vars[[N]]
    n <- prod(d)
    end <- i+n-1
    a <- array(i:end,d)
    i <- end+1
    ids[[N]] <- a
  }
  ids
### Named list of optimization variable ids, which will be used to
### construct the standard form matrices.
},ex=function(){
  vars <- make.ids(xi=5,beta0=1,beta=2)
})
constrain <- structure(function
### Make a constraint list. Write each constraint as c_1 v_1 + ... +
### c_n v_n (in)equality constant, where c_i are real coefficients and
### v_i are optimization variables.
(equality,
### "eq" or "geq" designating whether this is an equality or an
### inequality constraint.
 constant,
### constant on the right side of the inequality.
 ...
### Variables to include in this constraint, usually using something
### like variable("xi",i,j,-1).
 ){
  variables <- list(...)
  vars <- lapply(variables,structure,class=c("variable","list"))
  structure(list(variables=vars,equality=equality,constant=constant),
            class=c("constraint","list"))
},ex=function(){
  constraints <- list()
  for(i in 1:3){
    constraints <-
      c(constraints,list(constrain("geq",0,variable("xi",i))))
  }

  lapply(1:3,function(i)constrain("geq",0,variable("xi",i)))
})
variable <- structure(function
### Make a list that represents an optimization variable and its
### coefficient in the context of a constraint.
(name,
### Name of the optimization variable, created with make.ids.
 i=1,
### If the optimization variable is a vector, the element. If it is a
### matrix, the row.
 j=NULL,
### Column of the optimization variable.
 coef=1
### Coefficient of the optimization variable.
 ){
  structure(list(variable=name,i=i,j=j,coef=coef),
            class=c("variable","list"))
### List with elements "variable" "i" "j" "coef", of class "variable".
},ex=function(){
  variable("xi",2,coef=-1)
  variable("alpha",3,3)
})
print.variable <- function(x,...){
### print method that shows variable indices and coefficient
  index <- if(is.null(x$j)){
    sprintf("%d",x$i)
  }else{
    sprintf("%d,%d",x$i,x$j)
  }
  cat(sprintf("%10.3f*%s[%s]\n",x$coef,x$variable,index))
}
print.constraint <- function(x,...){
### print method that shows each variable and constant
  for(v in x$variables){
    print(v)
  }
  key <- c(geq=">=",eq="=")
  cat(key[x$equality],x$constant,"\n")
}
standard.form.constraints <- structure(function
### Convert constraints to standard form matrix and vector.
(constraints,
### List of constraints created using constrain().
 ids
### List of optimization variable ids created using make.ids().
 ){
  ## sort by type of constraint:
  eq <- sapply(constraints,"[[","equality")
  meq <- sum(eq=="eq") ## number of equality constraints
  ordered.constraints <- constraints[order(eq)]
  b0 <- sapply(ordered.constraints,"[[","constant")
  n.constraints <- length(ordered.constraints)
  n.variables <- length(unlist(ids))
  A <- matrix(0,n.variables,n.constraints)
  for(cn in seq_along(ordered.constraints)){
    constraint <- ordered.constraints[[cn]]
    for(vn in seq_along(constraint$variables)){
      v <- constraint$variables[[vn]]
      a <- ids[[v$variable]]
      if(is.null(a))stop(sprintf("variable %s not found",v$variable))
      A.row <- if(is.null(v$j)){
        a[v$i]
      }else{
        a[v$i,v$j]
      }
      A[A.row,cn] <- v$coef
    }
  }
  list(A=A,b0=b0,meq=meq)
### List of variable coefficients and constraint constants.
},ex=function(){
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
})
