# Newton's method

# f - function that calculates and returns the value of f at any given point x
# fprime - function that calculates and returns the derivative of f at any given point x
# fdoubleprime - function that calculates and returns the second derivative of f at any given point x
# x0 - initial starting point
# nIter - positive integer, number of iterations
NewtonMethod <- function(f, fprime, fdoubleprime, x0, nIter){
  
  # Implement Newton's method for nIter iterations, and return the vector of points xt, and function values f(xt)
  xvec <- rep(0, nIter)
  fvec <- rep(0, nIter)
  xvec[1] <- x0 - fprime(x0) / fdoubleprime(x0)
  fvec[1] < f(xvec[1])
  # At each iteration, save the function value f(x_t)
  for(i in 2:nIter) {
    xvec[i] <- xvec[i-1] - fprime(xvec[i-1]) / fdoubleprime(xvec[i-1])
    fvec[i] <- f(xvec[i])
  }
  return(list(xvec = xvec, fvec  = fvec))
}