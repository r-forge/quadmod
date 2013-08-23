library(quadmod)
vars <- make.ids(x=3)

## each variable must be positive.
positive <- vars$x >= 0
stopifnot(is.list(positive))
stopifnot(sapply(positive, is.constraint))
stopifnot(length(positive)==3)

## the sum of all x must be positive.
sum.pos <- sum(vars$x) >= 0
stopifnot(is.constraint(sum.pos))
stopifnot(length(sum.pos$variables)==3)

## the sum of all x must be 1.
sum.one <- sum(vars$x) == 1
stopifnot(is.constraint(sum.one))
stopifnot(length(sum.one$variables)==3)

## this linear combination must be equal to 1.
lc.one <- vars$x * c(-1,1,-2) == 1
stopifnot(is.constraint(lc.one))
stopifnot(length(lc.one$variables)==3)

## this variable must be positive.
one.pos <- vars$x[1] >= 0
stopifnot(is.constraint(one.pos))
stopifnot(length(one.pos$variables)==1)

## these two variables must be positive.
two.pos <- vars$x[c(1,3)] >= 0
stopifnot(sapply(two.pos, is.constraint))
stopifnot(length(two.pos)==2)

## each of these variables times these values must be positive.
each.lc.pos <- vars$x[] * c(-1,1,1) >= 0
stopifnot(is.list(each.lc.pos))
stopifnot(sapply(each.lc.pos, is.constraint))

## svm variable combination
svm <- make.ids(margin=1, weight=2)
svmc <- svm$weight*c(-1,1) + svm$margin*-1 >= 1
stopifnot(is.constraint(svmc))
stopifnot(length(svmc$variables)==3)

