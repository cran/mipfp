\name{mipfp-package}
\alias{mipfp-package}
\alias{mipfp}
\docType{package}
\title{
Multidimensional Iterative Proportional Fitting
}
\description{
An implementation of the iterative proportional
fitting (IPFP) procedure for updating an initial N-dimensional array (called a 
seed) with respect to given target marginal distributions. Those targets can 
also be multi-dimensional. This procedure is also able to estimate a 
(multi-dimensional) contingency table (encoded as an array) matching a given 
set of (multi-dimensional) margins. In that case, each cell of the seed must 
simply be set to 1.

The IPFP is also known as the RAS algorithm in economics and matrix raking or 
matrix scaling in computer science.
}
\details{
\tabular{ll}{
Package: \tab mipfp\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-26\cr
License: \tab GPL-2\cr
}
This package provides an implementation of the well-known iterative proportional
fitting procedure for updating a N-dimensional array with respect to given
target marginal distributions. Those targets can also multi-dimensional.

There is one function \code{Ipfp} documented in \link[mipfp]{Ipfp}.
}
\author{
Johan Barthelemy and Thomas Suesse with contributions from Mohammad Namazi-Rad.

Maintainer: Johan Barthelemy <johan@uow.edu.au>.
}
\references{
 Bacharach, M. (1965)
 Estimating Nonnegative Matrices from Marginal Data.
 \emph{International Economic Review} (Blackwell Publishing) 6 (3): 294-310.
 
 Deming, W. E., Stephan, F. F. (1940).
 On a Least Squares Adjustment of a Sampled Frequency Table When the Expected
 Marginal Totals are Known.
 \emph{Annals of Mathematical Statistics} 11 (4): 427-444.
 
 Stephan, F. F. (1942).
 Iterative method of adjusting frequency tables when expected margins are known.
 \emph{Annals of Mathematical Statistics} 13 (2): 166-178.
 
 Fienberg, S. E. (1970).
 An Iterative Procedure for Estimation in Contingency Tables.
 \emph{Annals of Mathematical Statistics} 41 (3): 907-917.
 
 Bishop, Y. M. M., Fienberg, S. E., Holland, P. W. (1975).
 \emph{Discrete Multivariate Analysis: Theory and Practice.}
 MIT Press. ISBN 978-0-262-02113-5.
}
\keyword{iterative proportional fitting procedure}
\keyword{ipfp}
\keyword{multiway contingency table update}
\keyword{RAS algorithm}
\keyword{matrix raking}
\keyword{matrix scaling}
\keyword{models}
\keyword{multivariate}
\keyword{array}

\seealso{
 \code{\link[ipfp]{ipfp}} for a package implementing ipfp to
 solve problems of the form Ax=b.
}
\examples{
# generation of an intial 2-way table to be updated
seed <- array(1,dim=c(2,2))
# desired targets (margins)
target.row <- c(87,13)
target.col <- c(52,48)
# storing the margins in a list
target.data <- list(target.col, target.row)
# storing the description of target data in a list
target.list <- list(1,2)
# calling the Ipfp function
res <- Ipfp(seed, target.list, target.data)
}
