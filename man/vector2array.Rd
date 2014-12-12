\name{vector2array}
\alias{vector2array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transforming a vector to an array
}
\description{
Transform a vector into a multidimensional array. The transformation is done
assuming that the last index of the array moves fastest. For instance, the
relation between a vector \eqn{v} of lenght 8 and an array \eqn{a} of dimensions
(2,2,2) is defined by \eqn{v = (a_{111}, a_{112}, a_{113}, a_{121}, a_{122},
\ldots,a_{333}}{ v = ( a[1,1,1], a[1,1,2], a[1,1,3], a[1,2,1], a[1,2,2],
\ldots,a[3,3,3] )}.
}
\usage{
vector2array(vect, dim.out)
}
\arguments{
  \item{vect}{
    The vector of length one or more to be transformed into an array.
}
  \item{dim.out}{
    The dimension attribute for the array to be created, that is an integer 
    vector of length one or more giving the maximal indices in each dimension.    
}
}
\value{
An array of dimensions given by \code{dim.out} filled with the data from the
input vector \code{vec}.
}
\author{
Thomas Suesse.

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\examples{
# generate a vector [1,2,...,27]
v <- seq(1:27)
# transform it into an array of dimension (3,3,3)
a <- vector2array(v,c(3,3,3))
}
\keyword{vector}
\keyword{array}
\keyword{manip}
