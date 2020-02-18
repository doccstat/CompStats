# Standardize X and Y
# X - n x p matrix of covariates
# Y - n x 1 response vector
standardizeXY <- function(X, Y){
  X <- as.matrix(X)
  Y <- as.vector(Y)
  n <- nrow(X)
  p <- ncol(X)
  if(n != length(Y)) { stop("Error: the row number of X are not equal to the length of Y!") }
  # Center Y
  Ymean <- mean(Y)
  Ytilde <- Y - Ymean

  # Center and scale X
  Xmeans <- colMeans(X)
  Xtilde <- X - matrix(Xmeans, n, p, byrow = T)
  fake_weights <- colSums(matrix(Xtilde^2, n, p))
  fake_weights[fake_weights == 0] <- n
  weights <- sqrt(n / fake_weights)
  Xtilde <- t(t(Xtilde) * weights)
  # Return the mean of Y and means of columns of X, as well as weights to be used in back-scaling (that is sqrt(X_j'X_j/n))
  return(list(Xtilde = Xtilde, Ytilde = Ytilde, Ymean = Ymean, Xmeans = Xmeans, weights = weights))
}

# Soft-thresholding of scalar a at level lambda
soft <- function(a, lambda){
  return(sign(a)*max(abs(a)-lambda,0))
}

# Calculates objective function of lasso given current values of Xtilde, Ytilde, beta and lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta - value of beta at which to evaluate the function
lasso <- function(Xtilde, Ytilde, beta, lambda){
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  # recompute n and p
  return(norm(Ytilde - Xtilde %*% beta, "F")^2 / (2*n) + lambda * sum(abs(beta)))
}

# Fit LASSO on standardized data for a given lambda
# Xtilde - centered and scaled X, n x p
# Ytilde - centered Y, n x 1
# lamdba - tuning parameter
# beta_start - p vector, an optiona starting point for coordinate-descent algorithm
# eps - precision level for convergence assessment, default 0.0001
fitLASSOstandardized <- function(Xtilde, Ytilde, lambda, beta_start = NULL, eps = 0.0001){
  # Check that n is the same between Xtilde and Ytilde
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  if(n != length(Ytilde)) { stop("Error: the row number of X are not equal to the length of Y!") }
  # Check that lambda is non-negative
  if(lambda < 0) { stop("Error: lambda should be non-negative!") }
  # Check for starting point beta_start. If none supplied, initialize with a vector of zeros. If supplied, check for compatibility with Xtilde in terms of p
  if(is.null(beta_start)) {
    beta_start <- as.vector(rep(0,p))
  } else {
    if(p != length(beta_start)) { 
      stop("Error: the length of beta_start is not equal to col num of Xtilde!") 
    }
  }
  # Coordinate-descent implementation. Stop when the difference between objective functions is less than eps.
  error <- 1
  fmin <- lasso(Xtilde, Ytilde, beta_start, lambda)
  residule <- Ytilde - Xtilde %*% beta_start

  while(error >= eps) {
    prev <- fmin
    for(i in 1:p) {
      temp <- soft(beta_start[i] + crossprod(Xtilde[ , i], residule)/n, lambda)
      residule <- residule + Xtilde[ , i] %*% (beta_start[i] - temp)
      beta_start[i] <- temp
    }
    fmin <- lasso(Xtilde, Ytilde, beta_start, lambda)
    error <- prev - fmin
  }
  beta <- beta_start
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
  n <- nrow(Xtilde)
  p <- ncol(Xtilde)
  if(n != length(Ytilde)) { stop("Error: the row number of X are not equal to the length of Y!") }
  # Check for the user-supplied lambda-seq (see below)
  # If lambda_seq is supplied, only keep values that are >=0, and make sure the values are sorted from largest to smallest. If none of the supplied values satisfy the requirement, print the warning message and proceed as if the values were not supplied.
  lambda_is_null <- is.null(lambda_seq)
  # flag whether we should create lambda_seq
  if(!lambda_is_null) {
    lambda_seq <- lambda_seq[lambda_seq >= 0]
    if(length(lambda_seq) > 0) {
      lambda_seq <- sort(lambda_seq, decreasing = TRUE)
      lambda_is_null <- FALSE
    }
    # else lambda_is_null remains TRUE 
  }
  # If lambda_seq is not supplied, calculate lambda_max (the minimal value of lambda that gives zero solution), and create a sequence of length n_lambda as
  if(lambda_is_null) {
    writeLines("Either lambda_seq is not provided or there are no valid values in it!\nSo The program will generate lambda_seq promptly.")
    lambda_max <- max(abs(crossprod(Xtilde, Ytilde))) / n
    if(lambda_max == 0) { lambda_max <- .Machine$double.xmin }
    lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))
  } else{
    n_lambda <- length(lambda_seq)
  }
  
  # Apply fitLASSOstandardized going from largest to smallest lambda (make sure supplied eps is carried over). Use warm starts strategy discussed in class for setting the starting values.
  beta_mat <- matrix(0, p, n_lambda)
  fmin_vec <- as.vector(rep(0, n_lambda))
  for(i in 1:n_lambda) {
    if(i == 1) { out <- fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[i], NULL , eps) }
    else { out <- fitLASSOstandardized(Xtilde, Ytilde, lambda_seq[i], beta_mat[ , i-1], eps) }
    beta_mat[ , i] <- out$beta
    fmin_vec[i] <- out$fmin
  }
  
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
  out <- standardizeXY(X, Y)
  Xtilde <- out$Xtilde
  Ytilde <- out$Ytilde
  Ymean <- out$Ymean
  Xmeans <- out$Xmeans
  weights <- out$weights

  writeLines("Warning from the author: I'm going the create two global variables ending in _Xingchi_Li+[0-9], which normally wouldn't overwrite your data, unless your name is Xingchi Li and your variables are the same as I defined.")
  assign("Xtilde_times_Ytilde_Xingchi_Li283482", crossprod(Xtilde, Ytilde), envir = .GlobalEnv)
  assign("Xtilde_times_Xtilde_Xingchi_Li234023", crossprod(Xtilde, Xtilde), envir = .GlobalEnv)
  # Fit Lasso on a sequence of values using fitLASSOstandardized_seq (make sure the parameters carry over)
  out <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq, n_lambda, eps)
  lambda_seq <- out$lambda_seq
  beta_mat <- out$beta_mat
  fmin_vec <- out$fmin_vec
  # Perform back scaling and centering to get original intercept and coefficient vector for each lambda

  beta_mat <- beta_mat * weights
  beta0_vec <- Ymean - Xmeans %*% beta_mat
  # Return output
  # lambda_seq - the actual sequence of tuning parameters used
  # beta_mat - p x length(lambda_seq) matrix of corresponding solutions at each lambda value (original data without center or scale)
  # beta0_vec - length(lambda_seq) vector of intercepts (original data without center or scale)
  return(list(lambda_seq = lambda_seq, beta_mat = beta_mat, beta0_vec = beta0_vec))
}

fitLASSO_mod2 <- function(X ,Y, lambda_seq = NULL, n_lambda = 50, eps = 0.0001){
  # Center and standardize X,Y based on standardizeXY function
  out <- standardizeXY(X, Y)
  Xtilde <- out$Xtilde
  Ytilde <- out$Ytilde
  Ymean <- out$Ymean
  Xmeans <- out$Xmeans
  weights <- out$weights

  # Fit Lasso on a sequence of values using fitLASSOstandardized_seq (make sure the parameters carry over)
  out <- fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq, n_lambda, eps)
  lambda_seq <- out$lambda_seq
  beta_mat <- out$beta_mat
  fmin_vec <- out$fmin_vec
  # Perform back scaling and centering to get original intercept and coefficient vector for each lambda

  beta_mat <- beta_mat * weights
  beta0_vec <- Ymean - Xmeans %*% beta_mat
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
  n <- nrow(X)
  p <- ncol(X)
  # Fit Lasso on original data using fitLASSO
  original_output <- fitLASSO(X, Y, lambda_seq, n_lambda, eps)
  lambda_seq <- original_output$lambda_seq
  beta_mat <- original_output$beta_mat
  beta0_vec <- original_output$beta0_vec
  n_lambda <- length(lambda_seq)

  # Split the data into K folds
  fold_identifier <- (sample(1:n) %% k) + 1
  # Calculate Lasso for each fold removed
  cv_lambda_seq_mat <- matrix(0, n_lambda, k)
  # n_lambda * k
  cv_beta_mat_mat <- matrix(0, p, n_lambda*k)
  # n_lambda * p*k
  cv_beta0_vec_mat <- matrix(0, n_lambda, k)
  # n_lambda * k

  for(i in 1:k) {
    # Calculate LASSO on that fold using fitLASSO
    Xtrain <- X[fold_identifier != i, ]
    Ytrain <- Y[fold_identifier != i]
    Xtest <- X[fold_identifier == i, ]
    Ytest <- Y[fold_identifier == i]
    out <- fitLASSO_mod2(Xtrain, Ytrain, lambda_seq, n_lambda, eps)
    # Any additional calculations that are needed for calculating CV and SE_CV(lambda)
    cv_lambda_seq_mat[ , i] <- out$lambda_seq
    cv_beta_mat_mat[ , ((i-1)*n_lambda+1):(i*n_lambda)] <- out$beta_mat
    cv_beta0_vec_mat[ , i] <- out$beta0_vec

  }
  # vectorization
  # cvm_first <- Y^2
  # cvm_second <- cv_beta0_vec_mat^2
  # cvm_third <- crossprod(X, cv_beta_mat_mat)

  cvm <- rep(0, n_lambda)
  # Calculate CV(lambda) and SE_CV(lambda) for each value of lambda
  cvm_element <- matrix(0, k, n_lambda)
  for(lambda_identifier in 1:n_lambda) {
    for(j in 1:k) {
      Ytest <- Y[fold_identifier == i]
      Xtest <- X[fold_identifier == i, ]
      cvm_element[j, lambda_identifier] <- sum((Ytest - cv_beta0_vec_mat[lambda_identifier, j] - Xtest %*% cv_beta_mat_mat[ , (j-1)*n_lambda+lambda_identifier]) ^2)
    }
  }
  
  # Find lambda_min
  lambda_min_identifier <- which.min(colSums(cvm_element))
  lambda_min <- lambda_seq[lambda_min_identifier]
  # Find lambda_1SE
  cvm <- colSums(cvm_element) / n
  for(j in 1:k) {
    cvm_element[j, ] = cvm_element[j, ] / sum(fold_identifier == j)
  }

  cvse <- rep(0, n_lambda)
  denominator <- sqrt(k)
  for(lambda_identifier in 1:n_lambda) {
    cvse[lambda_identifier] <- sd(cvm_element[ , lambda_identifier]) / denominator
  }
  temp <- cvm
  temp[temp > cvm[lambda_min_identifier] + cvse[lambda_min_identifier]] <- 0
  lambda_1se_identifier <- which.max(temp)
  lambda_1se <- lambda_seq[lambda_1se_identifier]
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

