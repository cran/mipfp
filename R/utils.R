# File mipfp/R/utils.R
# by Johan Barthelemy and Thomas Suesse
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# ------------------------------------------------------------------------------
# This file provides functions to convert array to vectors and conversely.
# ------------------------------------------------------------------------------

Array2Vector <- function(arr) {
  # Transform a N-dimensional array a to vector, where last index moves fastest.
  #
  # Author: T. Suesse
  #
  # Args:
  #   arr: The array to be transformed into a vector.
  #
  # Returns: A vector containing the input array data.
  
  # checking if an input array is specified
  if (is.null(arr) == TRUE) {
    stop('Error: no array arr specified!')
  }
  
  dim.array <- dim(arr)
  if (is.null(dim.array) == FALSE ) {
    arr <- aperm(arr, seq(length(dim.array), 1, by = -1))    
  }
  return(c(arr))
  
}

Vector2Array <- function(vect, dim.out) {
  # Transform a vector to an array with given dimensions, where last index moves 
  # fastest.
  #
  # Author: T. Suesse
  #
  # Args:
  #   vector: The vector to be transformed into an array.
  #   dim.out: The dimension of the generated array.
  #
  # Returns: An array filled with the input vector data.
  
  # checking if an input array is specified
  if (is.null(vect) == TRUE) {
    stop('Error: no vector vect specified!')
  }
  
  # checking if an input array is specified
  if (is.null(dim.out) == TRUE) {
    stop('Error: no dimension dim.out specified for the output array!')
  }
  
  l.dim.out<-length(dim.out)
  arr <- array(vect, dim.out[seq(l.dim.out,1,by=-1)])
  arr <- aperm(arr, seq(l.dim.out, 1, by = -1)) 
  return(arr)
  
}
