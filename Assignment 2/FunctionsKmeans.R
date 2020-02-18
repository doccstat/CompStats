# Function that implements K-means algorithm for a given number of iterations. The default number of iterations is 100.

MyKmeans <- function(X, K, M = NULL, numIter = 100) {
	# Check whether M is NULL or not. If NULL, initilize based on K random points from X. If not NULL, check for compatibility with X dimensions.
	X <- data.matrix(X)
	n <- nrow(X)
	if(is.null(M) || ncol(X) != ncol(M)) M <- as.matrix(X[sample.int(n, K), ])
	Y <- rep(0, n)
	# firstTerm <- rowSums(X^2)
	# firstTerm[i] == ||X[i, ]||^2
	while(numIter > 0) {
		numIter = numIter - 1
		Ypre <- Y
		secondTerm <- tcrossprod(X,M)
		# secondTerm[i,j] == ||X[i, ]||*||M[j, ]||
		thirdTerm <- rowSums(M^2)
		# thirdTerm[i] == ||M[i, ]||_2^2
		#for(i in 1:n) Y[i] <- which.min(firstTerm[i] - 2 * secondTerm[i, ] + thirdTerm)
		for(i in 1:n) Y[i] <- which.min(- 2 * secondTerm[i, ] + thirdTerm)
		# which.min tries to find the cluster
		if(identical(Y, Ypre)) numIter <- 0
		# check if Y converges
		if(ncol(X) == 1) {
			for(i in 1:K) M[i, ] <- mean(X[Y==i])
		} else {
			for(i in 1:K) M[i, ] <- colMeans(X[Y==i, ])
		} 
	}
	# Return the vector of assignments Y
	return(Y)
}