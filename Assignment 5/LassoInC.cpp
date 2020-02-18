#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double soft_c(double a, double lambda) {
	// Use ? : operator to speed the max operator
	return a < 0 ? -(-a-lambda > 0 ? -a-lambda : 0) : (a-lambda > 0 ? a-lambda : 0);
}

// [[Rcpp::export]]
double lasso_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& beta, double lambda) {
	int n = Xtilde.n_rows;
	return pow(norm(Ytilde - Xtilde * beta, "fro"), 2) / (2 * n) + lambda * norm(beta, 1);
}

// [[Rcpp::export]]
arma::colvec fitLASSOstandardized_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, double lambda, const arma::colvec& beta_start, double eps = 0.0001) {
	// Calculate n and p
	int n = Xtilde.n_rows, p = Xtilde.n_cols;
	// Predefine a big error
	double error = 1;
	double fmin = lasso_c(Xtilde, Ytilde, beta_start, lambda);
	arma::colvec residule = Ytilde - Xtilde * beta_start;
	// create beta since beta_start is const
	arma::colvec beta = beta_start;
	// Loop until error < eps
	while(error >= eps) {
		double prev = fmin;
		for(int i = 0; i < p; i++) {
			// save co-product to prevent redundant calculation
			double temp = soft_c(beta[i] + dot(Xtilde.col(i), residule) / n, lambda);
			residule = residule + Xtilde.col(i) * (beta[i] - temp);
			// update beta
			beta[i] = temp;
		}
		fmin = lasso_c(Xtilde, Ytilde, beta, lambda);
		// compute the fmin error
		error = prev - fmin;
	}
	return(beta);
}  

// [[Rcpp::export]]
arma::mat fitLASSOstandardized_seq_c(const arma::mat& Xtilde, const arma::colvec& Ytilde, const arma::colvec& lambda_seq, double eps = 0.0001) {
	int p = Xtilde.n_cols, n_lambda = lambda_seq.n_elem;
	arma::mat beta_mat(p, n_lambda);
	for(int i = 0; i < n_lambda; i++) {
		if(i == 0) { 
			// warm start need a starting point
			arma::colvec temp(p);
			temp.fill(0.0);
			beta_mat.col(i) = fitLASSOstandardized_c(Xtilde, Ytilde, lambda_seq[i], temp, eps); 
		}
		else { 
			// warm start
			beta_mat.col(i) = fitLASSOstandardized_c(Xtilde, Ytilde, lambda_seq[i], beta_mat.col(i-1), eps); 
		}
	}
	return(beta_mat);
}