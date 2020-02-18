# Assignment 2 - Implementation of the K-means clustering algorihtm

Please clone this directory locally, and use Rstudio preview for correct rendering of math equations.

## K-means

K-means is one of the most popular clustering algorithms. Given n data points $x_i$ in p dimensions, the algorithm iteratively divides the points into K clusters/groups such that the points within each group are most similar. Specifically, it aims to minimize
$$
\sum_{k=1}^K \sum_{x_i \in \text{cluster k}}\|x_i - \mu_k\|_2^2,
$$
where $\mu_k = \frac1{n_k}\sum_{x_i \in \text{cluster k}}x_i$ is the center or **centroid** of kth cluster with $n_k$ being the number of points in the cluster $k$. Here $\|x_i - \mu_k\|_2^2$ is the squared Euclidean distance between $x_i$ and $\mu_k$, although other distance metrics are possible.

## Algorithm

To minimize the above objective, K-means performs iterative adjustment of cluster centers. The algorithm is not guaranteed to find the absolute minimizer, but (hopefully) comes pretty close. To start, the algorithm chooses random $K$ points out of $n$ as the cluster centers $\mu_k$. At each iteration, the algorithm 
  
  * computes the distance from each point to each center $\mu_k$
  * assigns the points to the nearest center (among the K) 
  * after all $n$ points are assigned, the algorithm recomputes the centroid values $\mu_k$
  
The next iteration repeats these steps, and the process continues until some convergence criterion is met or the maximum number of iterations is reached.  

## Implementation

You are not allowed to use any external R libraries for this assignment, and you are not allowed to use __dist__ function.

You are provided the starter code in FunctionsKmeans.R that contains __MyKmeans(X, K, M, n_iter)__: This function takes n x p data matrix X, number of clusters K, (optional) initial K x p matrix of cluster centers M and (optional) pre-specified number of iterations for the algorithm, and returns the vector Y of length n of cluster assignments (numbers from 1 to K). More details are in the function comments. You are welcome to create additional functions if you want to.
  
## ZIPCODE example 
  
To test your function, we will use ZIPCODE.txt which provides 7,291 points of vectorized 16 by 16 pixels of image that should represent one of the 10 digits (from 0 to 9). The first column contains the correct cluster assignment (from 1 to 10), and the rest (256) are pixel values. The data is loaded in ZIPCODE_example.R, and you should use that R file for the analysis

Because the output of k-means algorithm is dependent on initial centroid values, you are asked to try the algorithm nRep = 50 times with different random initialization. Choose true number of clusters for K. At each replication, call your algorithm as implemented in __MyKmeans__ and evaluate the performance of the algorithm using the **RandIndex** implemented in the R package __fossil__. RandIndex takes values between 0 and 1, with 1 being the perfect match. The example code is provided but you will need to install the package __fossil__ first if you don't have it. At the end, report the mean RandIndex for your implementation across 50 replications.

## Grading for this assignment

Your assignment will be judged as follows

* correctness _(50% of the grade)_
* speed of implementation (you need to vectorize your code to pass the speed requirement) _(30% of the grade)_
* comments _(10% of the grade)_
* version control workflow (how readable are commits to others? was the project done in reasonable modules or in one big last minute push?) _(10% of the grade)_
* **(bonus)** you get 5 bonus points if your (correct) code is faster than my code
