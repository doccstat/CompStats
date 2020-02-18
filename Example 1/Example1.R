# One-dimensional optimization of a convex function

# Function definition and basic visualization
f <- function(x){
  (x - 50)^2 + exp(x)/50
}

plot(seq(-10,12, length = 100), f(seq(-10,12, length = 100)), type = "l", xlab = "x", ylab = "f(x)") 

# Use built in solvers
########################################################
# Using optimize to solve
optimize(f, interval = c(-100, 100))

# Function derivative
fprime <- function(x){
  2 * x - 100 + exp(x)/50
}

# Using uniroot to solve
uniroot(fprime, interval = c(-100, 100))

# Compare the two in terms of speed - optimize is faster here
library(microbenchmark)
microbenchmark(
  optimize(f, interval = c(-100, 100)),
  uniroot(fprime, interval = c(-100, 100))
)

# Can alternatively use optimize with derivative
gradient_difference <- function(x){
  return(fprime(x)^2)
}

optimize(gradient_difference,  interval = c(-100, 100))

microbenchmark(
  optimize(f, interval = c(-100, 100)),
  optimize(gradient_difference,  interval = c(-100, 100)),
  uniroot(fprime, interval = c(-100, 100))
)

# Use steepest descent with different step sizes
########################################################
nIter = 30
# Small step size
alpha = 0.001

out_small <- SteepestDescent(f, fprime, x0 = 0, alpha = alpha, nIter = nIter)

plot(1:nIter, out_small$xvec, type = 'o')

# Medium step size
alpha = 0.01

out_medium <- SteepestDescent(f, fprime, x0 = 0, alpha = alpha, nIter = nIter)

plot(1:nIter, out_medium$xvec, type = 'o')

# Large step size
alpha = 0.03

out_large <- SteepestDescent(f, fprime, x0 = 0, alpha = alpha, nIter = nIter)

plot(1:nIter, out_large$xvec, type = 'o')


# Use Newton method
########################################################
fdoubleprime <- function(x){
  2 + exp(x)/50
}

nIter = 30
# Starting point x0 = 0
out_Newton1 <- NewtonMethod(f, fprime, fdoubleprime, x0 = 0, nIter = nIter)
plot(1:nIter, out_Newton1$xvec, type = 'o')


# Starting point x0 = 5
out_Newton2 <- NewtonMethod(f, fprime, fdoubleprime, x0 = 5, nIter = nIter)
plot(1:nIter, out_Newton2$xvec, type = 'o')