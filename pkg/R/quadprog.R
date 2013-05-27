run.quadmod <- function
### Main function of the package. Should allow definition of the
### objective and constraints with a modeling language, then analyze
### and eliminate constraints and variables, then run solve.QP.
(vars,
### List of variables.
 objective,
### Some description of the objective function.
 constraints
### List of constraints.
 ){
  stop("not implemented, use run.quadprog for now")
}
  
run.quadprog <- structure(function
### Run quadprog with the specified constraints and annotate the
### output using the specified variables. The problem is: minimize_x
### -d'x + 1/2 x'Dx such that the constraints are verified.
(vars,
### List of variables.
 D,
### Quadratic term coefficients.
 d,
### Linear term coefficients.
 constraints
### List of constraints.
 ){
  info <- standard.form.constraints(constraints,vars)
  sol <- solve.QP(D,d,info$A,info$b0,info$meq)
  for(var.name in names(vars)){
    var.i <- vars[[var.name]]
    sol[[var.name]] <- sol$solution[var.i]
  }
  sol
},ex=function(){
  ## simplified linear svm example.
  set.seed(1)
  p <- 2
  y <- rep(c(-1,1),each=20)
  x <- replicate(p,rnorm(length(y),y))
  plot(x,col=y+2,asp=1)
  n <- nrow(x)

  ## Define the constraint functions.
  vars <- make.ids(slack=n,intercept=1,normal=p)
  constraints <- vars$slack >= 0
  for(i in 1:n){
    ivars <- with(vars,intercept*y[i] + sum(normal)*(x[i,]*y[i]) + slack[i])
    constraints <- c(constraints,list(ivars >= 1))
  }

  ## Define the objective function.
  n.vars <- length(unlist(vars))
  tolerance <- 1e-6
  Dvec <- rep(tolerance,n.vars)
  Dvec[vars$normal] <- 1
  D <- diag(Dvec)
  d <- rep(0,n.vars)
  d[vars$slack] <- -1 ## C == 1

  ## Convert to standard form and run the solver.
  sol <- run.quadprog(vars, D, d, constraints)

  ## Note that the optimal solution vectors can be accessed by their
  ## variable names.
  print(sol)
  stopifnot(all(names(vars) %in% names(sol)))

  with(sol,{
    title(paste("A linear Support Vector Machine (SVM):",
                "margin SVs circled, slack drawn in red for other SVs"))
    abline(-intercept/normal[2],-normal[1]/normal[2])
    abline((1-intercept)/normal[2],-normal[1]/normal[2],lty="dotted")
    abline((-1-intercept)/normal[2],-normal[1]/normal[2],lty="dotted")
    f <- function(x)intercept+sum(normal*x)
    yfx <- apply(x,1,f)*y
    on.margin <- abs(yfx-1)<tolerance
    points(x[on.margin,],cex=2,col=y[on.margin]+2)
    i <- yfx<1
    ## these 2 complicated formulas calculate the point on margin where
    ## the data point starts picking up slack
    x1 <- ((y[i]-intercept)*normal[1]-
           normal[1]*normal[2]*x[i,2]+
           normal[2]^2*x[i,1])/(normal[2]^2+normal[1]^2)
    x2 <- (y[i]-intercept)/normal[2]-normal[1]/normal[2]*x1
    segments(x[i,1],x[i,2],x1,x2,col="red")
  })
})

