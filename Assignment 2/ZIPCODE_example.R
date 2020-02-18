# Application of K-means algorithm to ZIPCODE data
source("FunctionsKmeans.R")
# Rand Index Calculation example
library(fossil)
cluster1 <- c(2,2,1,3)
cluster2 <- c(3,3,2,1)
rand.index(cluster1, cluster2) # clusters match

# Load the ZIPCODE data
zipcode <- read.table("ZIPCODE.txt", header = F)

# Extract the true digits
Y <- zipcode[ , 1]

# Extract the data points
X <- zipcode[ , -1]

# [ToDo] Try K-means algorithm nRep times with different starting points on ZIPCODE data. Calculate Rand Index at each replication
nRep <- 50

randIndex <- rep(0,nRep)

startTime <- Sys.time()

for(i in 1:nRep) {
	randIndex[i] <- rand.index(Y, MyKmeans(X = X, K = 10, M = NULL, numIter = 100))
}

endTime <- Sys.time()

# [ToDo] Report mean Rand Index

randIndex
mean(randIndex)
endTime - startTime