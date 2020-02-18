# Standardize X and Y
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X,Y){
  
  # Center Y
  
  # Center and scale X
  
  
  # Return the mean of Y and means of columns of X, as well as weights to be used in back-scaling (that is sqrt(X_j'X_j/n))
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# Soft-thresholding of vector a at level lambda
soft <- function(a,lambda){
  #Should return a vector of the same size as a
}

# Calculates objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
lasso <- function(Xtilde, Ytilde, beta, lambda){
 
}

# Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta_start - p vector, an optiona starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.0001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.0001){
  # Check that n is the same between Xtilde and Ytilde
  
  # Check that lambda is non-negative
  
  # Check for starting point beta_start. If none supplied, initialize with a vector of zeros. If supplied, check for compatibility with Xtilde in terms of p
  
  # Coordinate-descent implementation. Stop when the difference between objective functions is less than eps.
  
  # Return beta - the solution, fmin - optimal function value (value of objective at beta)
  return(list(beta = beta, fmin = fmin))
}

# Fit LASSO on standardized data for a sequence of lambda values. Sequence version of a previous function.
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.0001
fitLASSOstandardized_seq <- function(Xtilde, Ytilde, lambda_seq = NULL, n_lambda = 50, eps = 0.0001){
  # Check that n is the same between Xtilde and Ytilde
 
  # Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >=0, and make sure the values are sorted from largest to smallest. If none of the supplied values satisfy the requirement, print the warning message and proceed as if the values were not supplied.
  # If lambda_seq is not supplied, calculate lambda_max (the minimal value of lambda that gives zero solution), and create a sequence of length n_lambda as
  lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  
  # Apply fitLASSOstandardized going from largest to smallest lambda (make sure supplied eps is carried over). Use warm starts strategy discussed in class for setting the starting values.
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value
  # fmin_vec - length(lambda_seq) vector of corresponding objective function values at solution
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, fmin_vec = fmin_vec))
}

# Fit LASSO on original data using a sequence of lambda values
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# eps - precision level for convergence assessment, default 0.0001
fitLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 50, eps = 0.0001){
  # Center and standardize X,Y based on standardizeXY function
 
  # Fit Lasso on a sequence of values using fitLASSOstandardized_seq (make sure the parameters carry over)
 
  # Perform back scaling and centering to get original intercept and coefficient vector for each lambda
  
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}


# Fit LASSO and perform cross-validation to select the best fit
# X - n x p matrix of covariates
# Y - n x 1 response vector
# lamdba_seq - sequence of tuning parameters, optional
# n_lambda - length of desired tuning parameter sequence, is only used when the tuning sequence is not supplied by the user
# k - number of folds for k-fold cross-validation, default is 5
# eps - precision level for convergence assessment, default 0.0001
cvLASSO <- function(X ,Y, lambda_seq = NULL, n_lambda = 50, k = 5, eps = 0.0001){
  # Fit Lasso on original data using fitLASSO
 
  # Split the data into K folds
  
  # Calculate Lasso for each fold removed
  for (i in 1:k){
    # Calculate LASSO on that fold using fitLASSO
    
    # Any additional calculations that are needed for calculating CV and SE_CV(lambda)
  }
  # Calculate CV(lambda) and SE_CV(lambda) for each value of lambda
  
  # Find lambda_min

  # Find lambda_1SE
  
  # Return output
  # Output from fitLASSO on the whole data
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  # Other output
  # lambda_min - selected lambda based on minimal rule
  # lambda_1se - selected lambda based on 1SE rule
  # cvm - values of CV(lambda) for each lambda
  # cvse - values of SE_CV(lambda) for each lambda
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec, lambda_min = lambda_min, lambda_1se = lambda_1se, cvm = cvm, cvse = cvse))
}

