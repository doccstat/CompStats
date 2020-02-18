# Load the riboflavin data
require(hdi) # this should install hdi package if you don't have it already; otherwise library(hdi)
data(riboflavin) # this puts list with name riboflavin into the R environment, y - outcome, x - gene erpression
#dim(riboflavin$x) # n = 71 samples by p = 4088 predictors
#?riboflavin # this gives you more information on the dataset

# [ToDo] Source your lasso functions
source("LassoFunctions.R")
set.seed(2)
#X <- matrix(rnorm(100), 10, 10)
#Y <- matrix(rnorm(10), 10, 1)
#out1 <- standardizeXY(X, Y)
#all.equal(out$Xtilde, matrix(scale(X) * sqrt(nrow(X)/(nrow(X)-1)), nrow(X), ncol(X)))

#out2 <- fitLASSO(X, Y, 0, 50, 1e-4)
#outlm <- lm(Y ~ X)
#abs(c(out2$beta0_vec,out2$beta_mat)-as.vector(outlm$coefficients))

# Use your fitLASSO function on the riboflavin data with 30 tuning parameters

X <- riboflavin$x
Y <- riboflavin$y

out1 <- fitLASSO(X, Y, NULL, n_lambda = 30, 1e-4)

# Based on the above output, plot the number of non-zero elements in each beta versus the value of tuning parameter

beta_mat <- out1$beta_mat
yaxis <- colSums(beta_mat != 0)
xaxis <- out1$lambda_seq
plot(xaxis, yaxis)

# Use microbenchmark 10 times to check the timing of your fitLASSO function above with 30 tuning parameters

#library(microbenchmark)
#microbenchmark(fitLASSO(X, Y, NULL, n_lambda = 30, 1e-4), times = 10)

# Report your median timing in the comments here: (~5.5 seconds for Irina on her laptop)

# The median is 1.918103

# Use cvLASSO function on the riboflavin data with 30 tuning parameters

out2 <- cvLASSO(X, Y, NULL, n_lambda = 30, k = 5, eps = 1e-4)

# Based on the above output, plot the value of CV(lambda) versus tuning parameter

plot(xaxis, out2$cvm)



# test
source("LassoFunctions.R")
library(glmnet)

X1 <- matrix(c(-1, 3, 2, 1), 2, 2)
Y1 <- c(1, 0)

test11 <- fitLASSO(X1, Y1, NULL, n_lambda = 30, 1e-4)
test12 <- cvLASSO(X1, Y1, NULL, n_lambda = 30, k = 5, eps = 1e-4)

X2 <- matrix(c(-1, 3), 2, 1)
Y2 <- c(1, 0)

test21 <- fitLASSO(X2, Y2, 1/2, n_lambda = 1, 1e-4)
test22 <- cvLASSO(X2, Y2, NULL, n_lambda = 30, k = 5, eps = 1e-4)

X3 <- -1
Y3 <- 1

test31 <- fitLASSO(X3, Y3, NULL, n_lambda = 30, 1e-4)
test32 <- cvLASSO(X3, Y3, NULL, n_lambda = 30, k = 5, eps = 1e-4)

X4 <- matrix(c(-1, 2), 1, 2)
Y4 <- 1

test41 <- fitLASSO(X4, Y4, NULL, n_lambda = 30, 1e-4)
test42 <- cvLASSO(X4, Y4, NULL, n_lambda = 30, k = 5, eps = 1e-4)