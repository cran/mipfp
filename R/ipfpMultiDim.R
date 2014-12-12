# File mipfp/R/ipfpMultiDim.R
# by  Johan Barthelemy and Thomas Suesse
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

require(cmm)

Ipfp <- function(seed, target.list, target.data, print = FALSE, iter = 1000, 
                 tol = 1e-10, na.target = FALSE) {
  # Update an array using the iterative proportional fitting procedure.
  #
  # Author: J. Barthelemy
  #  
  # Args:
  #   seed: The initial multi-dimensional array to be updated. Each cell must
  #         be non-negative.
  #   target.list: A list of the target margins provided in target.data. Each
  #                component of the list is an array whose cells indicates
  #                 which dimension the corresponding margin relates to.
  #   target.data: A list containing the data of the target margins. Each
  #                component of the list is an array storing a margin.
  #                The list order must follow the one defined in target.list. 
  #                Note that the cells of the arrays must be non-negative, but
  #                can contains NA values.
  #   print: Verbose parameter: if TRUE prints the current iteration number
  #          and the value of the stopping criterion.
  #   iter: The maximum number of iteration allowed; must be greater than 0.
  #   tol: If the maximum absolute difference between two iteration is lower
  #        than the value specified by tol, then ipfp has reached convergence
  #        (stopping criterion); must be greater than 0.
  #   na.target: if set to TRUE, allows the targets to have NA cells. In that
  #              case the margins consistency is not checked.
  #
  # Returns: An array whose margins fit the target margins and of the same
  #          dimension as seed.
  
  # checking if NA in target cells if na.target is set to FALSE
  if (is.na(min(sapply(target.data, min))) & !na.target)  {
    stop('Error: NA values present in the margins - use na.target = TRUE!')
  }
  
  # checking non negativity condition for the seed and the target
  if (min(sapply(target.data, min), na.rm = na.target) < 0 | min(seed) < 0) {
    stop('Error: Target and Seed cells must be non-negative!')    
  }  
  
  # checking the strict positiviy of tol and iter
  if (iter < 1 | tol <= 0) {
    stop('Error: tol and iter must be strictly positive!')
  }
  
  # checking the margins consistency if no missing values in the targets
  check.margins <- TRUE
  
  if (na.target == FALSE) {
    if (length(target.data) > 1) {
      for (m in 2:length(target.data)) {      
        if (abs(sum(target.data[[m-1]]) - sum(target.data[[m]])) > 1e-10) {
          check.margins <- FALSE
          warning('Target not consistents - shifting to probabilities!
                  Check input data!\n')
          break
        }      
      }
    }
  } else {
    if (print) {
      cat('NOTE: Missing values present in target cells. ')
      cat('Margins consistency not checked!\n')  
    }        
  }
  
  # if margins are not consistent, shifting from frequencies to probabilities
  if (!check.margins) {
    seed <- seed / sum(seed)
    for (m in 1:length(target.data)) {
      target.data[[m]] <- target.data[[m]] / sum(target.data[[m]])
    }
  }
  
  if (print & check.margins & !na.target) {
    cat('Margins consistency checked!\n')
  } 
  
  # initial value is the seed
  result <- seed  
  converged <- FALSE
  tmp.evol.stp.crit <- vector(mode="numeric", length = iter)
  
  # ipfp iterations
  for (i in 1:iter) {
    
    if (print) {
      cat('... ITER', i, '\n')
    } 
    
    # saving previous iteration result (for testing convergence)
    result.temp <- result
            
    # loop over the constraints
    for (j in 1:length(target.list)) {
      # ... extracting current margins
      temp.sum <- apply(result, target.list[[j]], sum)
      # ... computation of the update factor, taking care of 0 and NA cells
      update.factor <- ifelse(target.data[[j]] == 0 | temp.sum == 0, 0,
                              target.data[[j]] / temp.sum)
      if (na.target == TRUE) {
        update.factor[is.na(update.factor)] <- 1;
      }
      # ... apply the update factor
      result <- sweep(result,target.list[[j]], update.factor, FUN = "*")
    }
    
    # stopping criterion
    stp.crit <- max(abs(result - result.temp))
    tmp.evol.stp.crit[i] <- stp.crit
    if (stp.crit < tol) {
      converged <- TRUE
      if (print) {
        cat('Convergence reached after', i, 'iterations!\n')
      } 
      break
    }
    
    if (print) {
      cat('       stoping criterion:', stp.crit, '\n')
    }
    
  }
  
  # checking the convergence
  if (converged == FALSE) {
    warning('IPFP did not converged after ', iter, ' iteration(s)! 
            This migh be due to 0 cells in the seed, maximum number 
            of iteration too low or tolerance too small\n')
  }        
  
  # computing final max difference between generated and target margins
  diff.margins <- vector(mode = "numeric", length = length(target.list))
  if (na.target == FALSE) {
    for (j in 1:length(target.list)) {
      diff.margins[j] = max(abs(target.data[[j]] 
                                - apply(result, target.list[[j]], sum))) 
    }
  }
  
  # storing the evolution of the stopping criterion
  evol.stp.crit <- tmp.evol.stp.crit[1:i]
  
  # gathering the results in a list
  result.list <- list("xi.hat" = result, "stp.crit" = stp.crit, 
                      "conv" = converged, "check.margins" = diff.margins,
                      "evol.stp.crit" = evol.stp.crit);
  
  # returning the result
  return(result.list)
  
}

array2vector<-function(arr) {
  # Transform a N-dimensional array a to vector, where last index moves fastest.
  #
  # Author: T. Suesse
  #
  # Args:
  #   arr: The array to be transformed into a vector.
  #
  # Returns: A vector containing the input array data.
  
  dim.array <- dim(arr)
  arr <- aperm(arr, seq(length(dim.array), 1, by = -1))
  
  return(c(arr))

}

vector2array<-function(vect, dim.out) {
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

  l.dim.out<-length(dim.out)
  arr <- array(vect, dim.out[seq(l.dim.out,1,by=-1)])
  arr <- aperm(arr, seq(l.dim.out, 1, by = -1))
  
  return(arr)
  
}

IpfpCovar <- function(estimate, seed, target.list, replace.zeros = 1e-10) {
  # Compute the covariance matrix of the estimators produced by Ipfp.
  #
  # This function determines the covariance matrix using the formula
  # given in the paper "Models for Contingency Tables With Known Margins When 
  # Target and Sampled Populations Differ" written by Little and Wu (1991).
  #
  # Author: J. Barthelemy
  #
  # Args:
  #   estimate: The array of estimate produced by the Ipfp function.
  #   seed: The initial multi-dimensional array updated by Ipfp. 
  #   target.list: A list of the target margins used by the Ipfp function. Each
  #                component of the list is an array whose cells indicates
  #                which dimension the corresponding margin relates to.
  #   replace.zeros: If 0-cells are to be found in either the seed or the
  #                  estimate arrays, then their values are replaced with this
  #                  value.
  #
  # Returns: A covariance matrix of the estimators (last index move fastest).
    
  n <- sum(seed)  
  seed.prob <- array2vector(seed / sum(seed))
  estimate.prob <- array2vector(estimate / sum(estimate))
  
  # checking if 0-cells values and replace them with a small value  
  seed.prob <- ifelse(seed.prob == 0, replace.zeros, seed.prob)    
  estimate.prob <- ifelse(estimate.prob == 0, replace.zeros, estimate.prob)
  
  # computation of the diagonal matrix filled with the inverses of seed and
  # estimated probabilities
  D.seed <- diag(1 / seed.prob)
  D.estimate <- diag(1 / estimate.prob)
  
  # computation of A such that A * vector(estimate) = vector(target.data)
  
  # ... one line filled with ones
  A.transp <- matrix(1, nrow = 1, ncol = length(estimate.prob))
  
  # ... constrainst (removing the first one since it is redundant information)
  for (j in 1:length(target.list)) {
    marg.mat <- cmm::MarginalMatrix(var = 1:length(dim(seed)), 
                                    marg = target.list[[j]], 
                                    dim = dim(seed))[-1,]
    A.transp <- rbind(marg.mat, A.transp, deparse.level = FALSE)
  }
  
  A <- t(A.transp)
  
  # computation of the orthogonal complement of A (using QR decomposition)
  K <- qr.Q(qr(A), complete = TRUE)[,(dim(A)[2]+1):dim(A)[1]]
  
  # computation of the variance  
  estimate.var <- (1 / n) * K %*% solve((t(K) %*% D.estimate %*% K)) %*%
                  t(K) %*% D.seed %*% K %*%
                  solve(t(K) %*% D.estimate %*% K) %*% t(K)
  
  # returning the result
  return(estimate.var)
  
}

