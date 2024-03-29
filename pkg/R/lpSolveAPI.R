run.lpSolveAPI <- structure(function
### Run lpSolve::lp to solve a linear program: maximize coefficients *
### variables such that the constraints are verified.
(variables,
### List of variables from make.ids.
 coefficients,
### Vector of objective function coefficients.
 constraints
### List of quadmod constraints.
 ){
  info <- standard.form.constraints(constraints, variables)
  n.vars <- nrow(info$A)  #info$A is variables x constraints
  stopifnot(length(coefficients) == n.vars)

  const.dir <- rep(">=", length(info$b0))
  if(info$meq > 0){
    const.dir[1:info$meq] <- "=="
  }

  ## Plug into lpSolveAPI.
  desc <- make.lp(ncol=n.vars)
  set.bounds(desc, lower=rep(-Inf, n.vars))
  for(i in 1:ncol(info$A)){
    add.constraint(desc, info$A[,i], const.dir[i], info$b0[i])
  }
  set.objfn(desc, -coefficients)
  status <- solve(desc)
  obj.const.vars <- get.primal.solution(desc)

  const.indices <- (1:length(constraints))+1
  solution <- obj.const.vars[(max(const.indices)+1):length(obj.const.vars)]
  fit <- list(objective=obj.const.vars[1],
              constraints=obj.const.vars[const.indices])
  for(var.name in names(variables)){
    var.indices <- variables[[var.name]]
    fit[[var.name]] <- solution[var.indices]
  }
  fit
},ex=function(){

  library(quadmod)
  data(good.bad.dist, envir=environment())
  feat <- good.bad.dist$features
  better <- good.bad.dist$better

  ## fit the linear max-margin comparison function to a set of
  ## pairs. good.bad.dist has features with 2 columns. The first
  ## column is x_i and the second column is x_i'. good.bad.dist has a
  ## vector better which gives y_i in {-1,0,1}, indicating which
  ## element of the pair is better: -1 means x_i is better, 1 means
  ## x_i' is better, and 0 means they are the same. The max margin
  ## comparison function is the solution to maximize_{mu,w} mu subject
  ## to mu < 1-|w(x_i'-x_i)|, for all i such that y_i=0, and for all
  ## other i, mu < -1 + w(x_i'-x_i)y_i. Translating this problem into
  ## standard form yields the following LP.

  vars <- make.ids(margin=1, weight=1)

  constraints <- list()
  for(i in 1:nrow(feat)){
    if(better[i] == 0){
      right.side <- -1
      yi.vec <- c(-1,1)
    }else{
      right.side <- 1
      yi.vec <- better[i]
    }
    for(yi in yi.vec){
      const <- with(vars,{
        weight*(feat[i,2]-feat[i,1])*yi + margin*-1 >= right.side
      })
      constraints <- c(constraints,list(const))
    }
  }

  n.vars <- length(unlist(vars))
  tolerance <- 1e-6
  Dvec <- rep(tolerance, n.vars)
  D <- diag(Dvec)
  d <- rep(0, n.vars)
  d[vars$margin] <- 1
  qp <- run.quadprog(vars, D, d, constraints)
  sol <- run.lpSolveAPI(vars, d, constraints)
  for(v.name in names(vars)){
    stopifnot(max(abs(qp[[v.name]]-sol[[v.name]])) < tolerance)
  }

  fxdiff <- sol$weight*(feat[,2]-feat[,1])
  thresh <- function(x)ifelse(x>1,1,ifelse(abs(x)<1,0,-1))
  ## check to make sure we have perfect prediction.
  stopifnot(thresh(fxdiff) == better)
  margin <- ifelse(better==0,{
    1-abs(fxdiff)
  },{
    -1 + better * fxdiff
  })
  on.margin <- abs(margin - sol$margin)<tolerance
  margin.points <- feat[on.margin,]
  margin.better <- better[on.margin]
  boundary <- ifelse(margin.better==0,{
    ifelse(margin.points[,2]>margin.points[,1], -1, 1)
  },{
    ifelse(margin.better == 1, -1, 1)
  })

  boundary.x <- margin.points[,2]
  boundary.y <- boundary.x + boundary/sol$weight
  margin.df <- data.frame(margin.points,boundary.x,boundary.y)
  point.df <- data.frame(feat,better=factor(better))
  line.df <- data.frame(slope=1,intercept=c(-1,1)/sol$weight)

  library(ggplot2)
  p <- ggplot(,aes(X2,X1))+
    geom_point(aes(colour=better), data=point.df)+
  coord_equal()+
  geom_abline(aes(slope=slope,intercept=intercept),data=line.df)+
  scale_colour_discrete("$y_i$")+
  xlab("$x_i'$")+
  ylab("$x_i$")+
  geom_segment(aes(xend=boundary.x, yend=boundary.y), data=margin.df)

  print(p)
})
