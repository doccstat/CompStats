# Steepest descent with the fixed number of iterations

# f - function that calculates and returns the value of f at any given point x
# fprime - function that calculates and returns the derivative of f at any given point x
# x0 - initial starting point
# alpha - positive number, step size
# nIter - positive integer, number of iterations
SteepestDescent <- function(f, fprime, x0, alpha, nIter){
  
  # Perform steepest descent update for nIter iterations
  xvec <- rep(0, nIter)
  fvec <- rep(0, nIter)
  xvec[1] <- x0 - alpha * fprime(x0)
  fvec[1] < f(xvec[1])
  # At each iteration, save the function value f(x_t)
  for(i in 2:nIter) {
    xvec[i] <- xvec[i-1] - alpha * fprime(xvec[i-1])
    fvec[i] <- f(xvec[i])
  }
  # Return the vector of x values, as well as the vector of function values across iterations
  return(list(xvec = xvec, fvec  = fvec))
}


