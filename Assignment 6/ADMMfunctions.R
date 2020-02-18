# Matrix soft-thresholding for ell_1 norm
# A - n x p matrix
# lambda - non-negative parameter to threshold
soft <- function(A, lambda){
  return(sign(A) * ((abs(A)-lambda) * (abs(A)-lambda > 0)))
  # Should return the matrix whose elements were soft-thresholded at level lambda
}

# Matrix singular values soft-thresholding for nuclear norm
# A - n x p matrix
# lambda - non-negative parameter to threshold
soft_nuclear <- function(A, lambda){
  svd_A <- svd(A)
  # Fill in
  newd <- svd_A$d

  newA <- tcrossprod(svd_A$u %*% soft(diag(newd), lambda), svd_A$v)
  # Should return a list with
  # newd - soft-thresholded singular values, vector
  # newA - n x p matrix with singular values newd
  return(list(newd = newd, newA = newA))
}

# Objective value of robustPCA, |L|_* + gamma|S|_1,
# Ld - vector of singular values of L (this avoids recalculating svd more than needed)
# S - current value of matrix S
# gamma - current value of gamma
robustPCAobj <- function(Ld, S, gamma = 0.1){
  # Fill in, should return the scalar - value of the objective
  return(sum(svd(Ld)$d) + gamma * sum(abs(S)))
}


# ADMM algorithm for Robust PCA
# M - given matrix, n by p
# gamma - positive parameter, default value 0.1
# Sinit - starting value for S, if none is supplied initialize with matrix of zeros. If supplied, check for compatibility
# etainit - starting value of eta, if none is supplied initialize with zeros. If supplied, check for compatibility.
# tau - ADMM parameter (scaled ADMM version), default value is 1
# eps - convergence tolerance criteria (difference in objective function values), default value is 0.0001
robustPCAadmm <- function(M, gamma = 0.1, Sinit = NULL, etainit = NULL, tau = 1, eps = 0.0001){
  # Initialize S and eta, and do necessary compatibility checks
  n <- nrow(M)
  p <- ncol(M)
  
  if(is.null(Sinit)) {
    Sinit <- matrix(rep(0, n*p), n, p)
  } else {
    if(n != nrow(Sinit)) { stop("row number of Sinit doesn't match") }
    if(p != ncol(Sinit)) { stop("col number of Sinit doesn't match") }
  }
  if(is.null(etainit)) {
    etainit <- matrix(rep(0, n*p), n, p)
  } else {
    if(n != nrow(etainit)) { stop("row number of etainit doesn't match") }
    if(p != ncol(etainit)) { stop("col number of etainit doesn't match") }
  }
  S <- Sinit
  eta <- etainit
  # Initialize L as L = M - S
  L <- M - S
  ## Implement ADMM algorithm, which alternates updates of L, S and eta until convergence
  error <- 1
  prev <- robustPCAobj(L, S, gamma)
  while(error > eps) {
    L <- soft_nuclear(M - S - eta, tau)$newA
    S <- soft(M - L - eta, gamma * tau)
    eta <- eta + L + S - M
    obj <- robustPCAobj(L, S, gamma)
    error <- prev - obj
    prev <- obj
  }
  # Return L, S and eta at convergence
  return(list(L = L, S = S, eta = eta))
}