# Function that implements multi-class logistic regression.
#############################################################
# Description of supplied parameters:
# X - n x p training data, 1st column should be 1s to account for intercept
# Y - a vector of size n of class labels, from 0 to K-1
# Xt - ntest x p testing data, 1st column should be 1s to account for intercept
# Yt - a vector of size ntest of test class labels, from 0 to K-1
# numIter - number of FIXED iterations of the algorithm, default value is 50
# eta - learning rate, default value is 0.1
# lambda - ridge parameter, default value is 0.1
# beta_init - (optional) initial starting values of beta for the algorithm, should be p x K matrix 

## Return output
##########################################################################
# beta - p x K matrix of estimated beta values after numIter iterations
# error_train - (numIter + 1) length vector of training error % at each iteration (+ starting value)
# error_test - (numIter + 1) length vector of testing error % at each iteration (+ starting value)
# objective - (numIter + 1) length vector of objective values of the function that we are minimizing at each iteration (+ starting value)
LRMultiClass <- function(X, Y, Xt, Yt, numIter = 50, eta = 0.1, lambda = 0.1, beta_init = NULL) {
  ## Check the supplied parameters
  ###################################
  # Check that the first column of X and Xt are 1s, if not - display appropriate message and stop execution.
  X <- matrix(as.matrix(X),ncol=ncol(X),dimnames=NULL)
  n <- length(Y)
  p <- ncol(X)
  ntest <- length(Yt)
  if(!identical(X[ , 1], as.vector(rep(1,n)))) { stop("The first column should be 1s!") }
  # Check for compatibility between X and Y
  # Check for compatibility between Xt and Yt
  # Check for compatibility between X and Xt
  stopifnot(nrow(X) == n && p == ncol(Xt) && nrow(Xt) == ntest)
  # Check eta is positive
  # Check lambda is non-negative
  stopifnot(eta > 0 && lambda >= 0)
  # Check whether beta_init is NULL. If NULL, initialize beta with p x K matrix of zeroes and construct corresponding pbeta. If not NULL, check for compatibility with what has been already supplied.
  K <- max(Y) + 1
  if(K < max(Yt) + 1) {
    stop("labels of Y should not be less than Yt!")
  }
  if(is.null(beta_init)) { 
    beta_init <- matrix(rep(0,p*K),p,K) 
  } else { 
    K <- max(K, ncol(beta_init))
    stopifnot(nrow(beta_init) == p) 
  }
  ## Calculate corresponding pk, objective value at the starting point, training error and testing error given the starting point
  ##########################################################################
  # Compute the objective value, training error and testing error in the for loop is better since in R the index starts from 1
  
  ## Newton's method cycle - implement the update EXACTLY numIter iterations
  ##########################################################################
  # Within one iteration: perform the update, calculate updated objective function and training/testing errors

  # Declaring the variables to avoid reallocation memory
  error_train <- rep(0, numIter+1)
  error_test <- error_train
  objective <- error_train

  # Declare matrices in advance for better performance

  for(iter in 1:(numIter+1)) {
    exponent <- exp(X %*% beta_init)
    rowsums_exp <- rowSums(exponent)
    P <- exponent / rowsums_exp
    W <- P - P^2

    # Construct matrices P and W

    label_train <- apply(X %*% beta_init, 1, which.max) - 1
    label_test <- apply(Xt %*% beta_init, 1, which.max) - 1

    # Predict labels based on current beta

    objective[iter] <- - sum(colSums(t(X)*beta_init[ , Y+1])) + sum(log(rowsums_exp)) + lambda*norm(beta_init,"F")/2
    error_train[iter] <- sum(label_train != Y) / n
    error_test[iter] <- sum(label_test != Yt) / ntest
    for(k in 1:K) {

      # Construct matrix of identity

      temp <- matrix(rep(0,n), n, 1)
      temp[Y == (k-1)] <- 1

      # Update beta

      beta_init[ , k] <- beta_init[ , k] - eta * solve(crossprod(X*W[ , k],X)+diag(rep(lambda,p)),crossprod(X,P[ , k]-temp)+lambda*beta_init[ , k])
    }
  }
  ## Return output
  ##########################################################################
  # beta - p x K matrix of estimated beta values after numIter iterations
  # error_train - (numIter + 1) length vector of training error % at each iteration (+ starting value)
  # error_test - (numIter + 1) length vector of testing error % at each iteration (+ starting value)
  # objective - (numIter + 1) length vector of objective values of the function that we are minimizing at each iteration (+ starting value)
  return(list(beta = beta_init, error_train = error_train, error_test = error_test, objective =  objective))
}