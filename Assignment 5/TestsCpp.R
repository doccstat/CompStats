# Change this to source your lasso functions
source("LassoFunctions.R")

library(microbenchmark)
library(hdi)
data(riboflavin)

# Header for Rcpp and RcppArmadillo
library(Rcpp)
library(RcppArmadillo)

# Source your C++ funcitons
sourceCpp("LassoInC.cpp")

X <- riboflavin$x
Y <- riboflavin$y
out <- standardizeXY(X, Y)
Xtilde <- out$Xtilde
Ytilde <- out$Ytilde

X1 <- matrix(c(-1, 3, 2, 1), 2, 2)
Y1 <- c(1, 0)
out <- standardizeXY(X1, Y1)
Xtilde1 <- out$Xtilde
Ytilde1 <- out$Ytilde


X2 <- matrix(c(-1, 3), 2, 1)
Y2 <- c(1, 0)
out <- standardizeXY(X2, Y2)
Xtilde2 <- out$Xtilde
Ytilde2 <- out$Ytilde

X3 <- -1
Y3 <- 1
out <- standardizeXY(X3, Y3)
Xtilde3 <- out$Xtilde
Ytilde3 <- out$Ytilde

X4 <- matrix(c(-1, 2), 1, 2)
Y4 <- 1
out <- standardizeXY(X4, Y4)
Xtilde4 <- out$Xtilde
Ytilde4 <- out$Ytilde

n_lambda <- 30

lambda_max <- max(abs(crossprod(Xtilde, Ytilde))) / nrow(Xtilde)
if(lambda_max == 0) { lambda_max <- .Machine$double.xmin }
lambda_seq = exp(seq(log(lambda_max), log(0.01), length = n_lambda))

lambda_max <- max(abs(crossprod(Xtilde1, Ytilde1))) / nrow(Xtilde1)
if(lambda_max == 0) { lambda_max <- .Machine$double.xmin }
lambda_seq1 = exp(seq(log(lambda_max), log(0.01), length = n_lambda))

lambda_max <- max(abs(crossprod(Xtilde2, Ytilde2))) / nrow(Xtilde2)
if(lambda_max == 0) { lambda_max <- .Machine$double.xmin }
lambda_seq2 = exp(seq(log(lambda_max), log(0.01), length = n_lambda))

lambda_max <- max(abs(crossprod(Xtilde3, Ytilde3))) / nrow(Xtilde3)
if(lambda_max == 0) { lambda_max <- .Machine$double.xmin }
lambda_seq3 = exp(seq(log(lambda_max), log(0.01), length = n_lambda))

lambda_max <- max(abs(crossprod(Xtilde4, Ytilde4))) / nrow(Xtilde4)
if(lambda_max == 0) { lambda_max <- .Machine$double.xmin }
lambda_seq4 = exp(seq(log(lambda_max), log(0.01), length = n_lambda))

beta_start <- as.vector(rep(0,ncol(Xtilde)))
beta_start1 <- as.vector(rep(0,ncol(Xtilde1)))
beta_start2 <- as.vector(rep(0,ncol(Xtilde2)))
beta_start3 <- as.vector(rep(0,ncol(Xtilde3)))
beta_start4 <- as.vector(rep(0,ncol(Xtilde4)))

# Do 2 tests for soft-thresholding function below
#################################################

all.equal(soft_c(1,-3), soft(1,-3))
all.equal(soft_c(1, 0), soft(1, 0))
all.equal(soft_c(-2,-1), soft(-2,-1))
all.equal(soft_c(-2, 3), soft(-2, 3))

# Do 2 tests for lasso objective function below
#################################################

all.equal(lasso(Xtilde, Ytilde, beta_start, lambda = lambda_seq[n_lambda/2]), lasso_c(Xtilde, Ytilde, beta_start, lambda = lambda_seq[n_lambda/2]))
all.equal(lasso(Xtilde1, Ytilde1, beta_start1, lambda = lambda_seq1[n_lambda/2]), lasso_c(Xtilde1, Ytilde1, beta_start1, lambda = lambda_seq1[n_lambda/2]))
all.equal(lasso(Xtilde2, Ytilde2, beta_start2, lambda = lambda_seq2[n_lambda/2]), lasso_c(Xtilde2, Ytilde2, beta_start2, lambda = lambda_seq2[n_lambda/2]))
all.equal(lasso(Xtilde3, Ytilde3, beta_start3, lambda = lambda_seq3[n_lambda/2]), lasso_c(Xtilde3, Ytilde3, beta_start3, lambda = lambda_seq3[n_lambda/2]))
all.equal(lasso(Xtilde4, Ytilde4, beta_start4, lambda = lambda_seq4[n_lambda/2]), lasso_c(Xtilde4, Ytilde4, beta_start4, lambda = lambda_seq4[n_lambda/2]))

# Do 2 tests for fitLASSOstandardized function below
#################################################

all.equal(fitLASSOstandardized(Xtilde, Ytilde, lambda = lambda_seq[n_lambda/2], beta_start, eps = 0.0001)$beta, fitLASSOstandardized_c(Xtilde, Ytilde, lambda = lambda_seq[n_lambda/2], beta_start, eps = 0.0001))
all.equal(fitLASSOstandardized(Xtilde1, Ytilde1, lambda = lambda_seq1[n_lambda/2], beta_start1, eps = 0.0001)$beta, fitLASSOstandardized_c(Xtilde1, Ytilde1, lambda = lambda_seq1[n_lambda/2], beta_start1, eps = 0.0001))
all.equal(fitLASSOstandardized(Xtilde2, Ytilde2, lambda = lambda_seq2[n_lambda/2], beta_start2, eps = 0.0001)$beta, fitLASSOstandardized_c(Xtilde2, Ytilde2, lambda = lambda_seq2[n_lambda/2], beta_start2, eps = 0.0001))
all.equal(fitLASSOstandardized(Xtilde3, Ytilde3, lambda = lambda_seq3[n_lambda/2], beta_start3, eps = 0.0001)$beta, fitLASSOstandardized_c(Xtilde3, Ytilde3, lambda = lambda_seq3[n_lambda/2], beta_start3, eps = 0.0001))
all.equal(fitLASSOstandardized(Xtilde4, Ytilde4, lambda = lambda_seq4[n_lambda/2], beta_start4, eps = 0.0001)$beta, fitLASSOstandardized_c(Xtilde4, Ytilde4, lambda = lambda_seq4[n_lambda/2], beta_start4, eps = 0.0001))

# Do microbenchmark on fitLASSOstandardized vs fitLASSOstandardized_c
######################################################################

microbenchmark(
	fitLASSOstandardized(Xtilde, Ytilde, lambda = lambda_seq[n_lambda/2], beta_start, eps = 0.0001),
	fitLASSOstandardized_c(Xtilde, Ytilde, lambda = lambda_seq[n_lambda/2], beta_start, eps = 0.0001),
	times = 10
)

 #   fitLASSOstandardized(Xtilde, Ytilde, lambda = lambda_seq[n_lambda/2],      beta_start, eps = 1e-04) 6340.0564
 # fitLASSOstandardized_c(Xtilde, Ytilde, lambda = lambda_seq[n_lambda/2],      beta_start, eps = 1e-04)   25.4199
 #        lq       mean    median        uq       max neval
 # 6351.2956 6437.02940 6395.4371 6467.8794 6776.8111    10
 #   26.2544   28.17071   26.6652   28.6146   39.3629    10

# Do 2 tests for fitLASSOstandardized_seq function below
#################################################

all.equal(fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq, n_lambda = 30, eps = 0.0001)$beta_mat, fitLASSOstandardized_seq_c(Xtilde, Ytilde, lambda_seq = lambda_seq, eps = 0.0001))
all.equal(fitLASSOstandardized_seq(Xtilde1, Ytilde1, lambda_seq = lambda_seq1, n_lambda = 30, eps = 0.0001)$beta_mat, fitLASSOstandardized_seq_c(Xtilde1, Ytilde1, lambda_seq = lambda_seq1, eps = 0.0001))
all.equal(fitLASSOstandardized_seq(Xtilde2, Ytilde2, lambda_seq = lambda_seq2, n_lambda = 30, eps = 0.0001)$beta_mat, fitLASSOstandardized_seq_c(Xtilde2, Ytilde2, lambda_seq = lambda_seq2, eps = 0.0001))
all.equal(fitLASSOstandardized_seq(Xtilde3, Ytilde3, lambda_seq = lambda_seq3, n_lambda = 30, eps = 0.0001)$beta_mat, fitLASSOstandardized_seq_c(Xtilde3, Ytilde3, lambda_seq = lambda_seq3, eps = 0.0001))
all.equal(fitLASSOstandardized_seq(Xtilde4, Ytilde4, lambda_seq = lambda_seq4, n_lambda = 30, eps = 0.0001)$beta_mat, fitLASSOstandardized_seq_c(Xtilde4, Ytilde4, lambda_seq = lambda_seq4, eps = 0.0001))

# Do microbenchmark on fitLASSOstandardized_seq vs fitLASSOstandardized_seq_c
######################################################################

microbenchmark(
	fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq, n_lambda = 30, eps = 0.0001),
	fitLASSOstandardized_seq_c(Xtilde, Ytilde, lambda_seq = lambda_seq, eps = 0.0001),
	times = 10
)

 # fitLASSOstandardized_seq(Xtilde, Ytilde, lambda_seq = lambda_seq,      n_lambda = 30, eps = 1e-04) 15889.6283
 #              fitLASSOstandardized_seq_c(Xtilde, Ytilde, lambda_seq = lambda_seq,      eps = 1e-04)    75.0632
 #         lq        mean      median         uq        max neval
 # 16049.0580 16449.36591 16286.09450 16649.8337 17647.5987    10
 #    76.8499    80.08555    77.90205    80.1547    98.8556    10

# Tests on riboflavin data
##########################
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression

# Separate standardization from lasso sequence
out <- standardizeXY(riboflavin$x, riboflavin$y)
# This is just to create lambda_seq, can be done faster, but this is simpler for now
outl <- fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, n_lambda = 30)

# The code below should assess your speed improvement on riboflavin data
microbenchmark(
	fitLASSOstandardized_seq(out$Xtilde, out$Ytilde, outl$lambda_seq),
	fitLASSOstandardized_seq_c(out$Xtilde, out$Ytilde, outl$lambda_seq),
	times = 10
)
