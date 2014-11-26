Ipfp <- function(seed, target.list, target.data, print = FALSE, iter = 100, 
                 tol = 1e-10) {
  # Update an array using the iterative proportional fitting procedure.
  #  
  # Args:
  #   seed: The initial multi-dimensional array to be updated. Each cell must
  #         be greater than 0.
  #   target.list: A list of the target margins provided in target.data. Each
  #                component of the list is an array whose cells indicates which
  #                dimension the corresponding margin relates to.
  #   target.data: A list containing the data of the target margins. Each
  #                component of the list is an array storing a margin.
  #                The list order must follow the one defined in target.list. 
  #                Note that the cells of the arrays must be greater than 0.
  #   print: Verbose parameter: if TRUE prints the current iteration number
  #          and the value of the stopping criterion.
  #   iter: The maximum number of iteration allowed; must be greater than 0.
  #   tol: If the maximum absolute difference between two iteration is lower
  #        than the value specified by tol, then ipfp has reached convergence
  #        (stopping criterion); must be greater than 0.
  #
  # Returns: An array whose margins fit the target margins and of the same
  #          dimension as seed.
  
  # checking non negativity condition for the seed and the target
  if (min(sapply(target.data, min)) < 0 | min(seed) < 0) {
    stop('Error: Target and Seed cells must be non-negative!')    
  }
    
  # checking the strict positiviy of tol and iter
  if (iter < 1 | tol <= 0) {
    stop('Error: tol and iter must be strictly positive!')
  }
  
  # checking the margins consistency
  check.margins <- TRUE
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
  
  # if margins are not consistent, shifting from frequencies to probabilities
  if (!check.margins) {
    seed <- seed / sum(seed)
    for (m in 1:length(target.data)) {
      target.data[[m]] <- target.data[[m]] / sum(target.data[[m]])
    }
  }
  
  if (print & check.margins) {
    cat('Margins consistency checked\n')
  } 
  
  # initial value is the seed
  result <- seed  
  converged <- FALSE
  
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
      # ... computation of the update factor, taking care of 0 cells
      update.factor <- ifelse(target.data[[j]] == 0 | temp.sum == 0, 0,
                              target.data[[j]] / temp.sum)
      # ... apply the update factor
      result <- sweep(result,target.list[[j]], update.factor, FUN = "*")
    }
    
    # stopping criterion
    stp.crit <- max(abs(result - result.temp))
    if (stp.crit < tol) {
      converged <- TRUE
      if (print) {
        cat('Convergence reached after', i, 'iterations!\n')
      } 
      break
    }
    
    if (print) {
      cat ('       stoping criterion:', stp.crit, '\n')
    }
    
  }
  
  # checking the convergence
  if (converged == FALSE) {
    warning('IPFP did not converged after ', iter, ' iteration(s)! 
            This migh be due to 0 cells in the seed, maximum number 
            of iteration too low or tolerance too small\n')
  }  
  
  # returning the result
  return(result)
  
}
