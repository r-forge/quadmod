\name{good.bad.dist}
\alias{good.bad.dist}
\docType{data}
\title{
Labeled pairs of features
}
\description{
A set of 20 layouts of a graph with 2 nodes was created. For every
layout we calculate the Euclidean distance between the two nodes. All
pairs were considered, and the distances of a pair of layouts appears as
a row of the features matrix. The better vector gives a label for each
pair: -1 if the left layout (column 1) is better, 1 if the right layout
(column 2) is better, and 0 if they are the same. The goal is to learn a
function that gives a score to each layout that respects the labels.
}
\usage{data(good.bad.dist)}
\format{
  The format is:
List of 2
 $ features: num [1:190, 1:2] 202 202 202 202 202 ...
 $ better  : num [1:190] 0 -1 -1 0 -1 -1 0 -1 0 0 ...
}
\references{
Learning to compare graph layouts, TD Hocking and S Spanurattana.
}
\keyword{datasets}
