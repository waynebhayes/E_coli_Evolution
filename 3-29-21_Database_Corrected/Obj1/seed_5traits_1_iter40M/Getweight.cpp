#include <fstream>
#include <iostream>
#include <algorithm>
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
using namespace std;
using namespace Rcpp;
using namespace RcppEigen;

//[[Rcpp::export]]
NumericMatrix MatrixMul(NumericMatrix x, NumericMatrix y);
NumericMatrix MatrixMul_schur(NumericMatrix x, NumericMatrix y);  
NumericVector Getweight(NumericVector GR_tr);

//[[Rcpp::export]]
NumericVector Getweight(NumericVector GR_tr){
	cout<<"BEGIN"<<endl;
	NumericVector GR_tr_sorted = clone(GR_tr).sort();
	NumericVector GR_tr_count (unique(GR_tr_sorted).length(), 0);
	//Rcout<<as<arma::rowvec>(GR_tr[seq_len(20)])<<endl;
	//cout<<GR_tr_sorted<<endl;
	NumericVector GR_tr_weights (GR_tr.length(), 0);
	float weight;
	int pos = 0;
	double current = GR_tr_sorted[0];
	
	for (int i=0; i<GR_tr_sorted.length(); i++){
		if (GR_tr_sorted[i] == current){
			GR_tr_count[pos]++;
		} else {
			current = GR_tr_sorted[i];
			pos++;
			GR_tr_count[pos]++;
		}	
	}
	GR_tr_sorted = unique(GR_tr_sorted).sort();
	for (int i=0; i<GR_tr.length(); i++ ) {
		pos = distance(GR_tr_sorted.begin(), std::find(GR_tr_sorted.begin(), GR_tr_sorted.end(), GR_tr[i]));
		if (pos == 0){
			weight = (GR_tr_sorted[1] - GR_tr_sorted[0]) / (GR_tr_count[0] + GR_tr_count[1] - 1);
		}
		else if (pos == GR_tr_sorted.length() - 1){
			weight = (GR_tr_sorted[GR_tr_sorted.length() - 1] - GR_tr_sorted[GR_tr_sorted.length() - 2]) / (GR_tr_count[GR_tr_sorted.length() - 1] + GR_tr_count[GR_tr_sorted.length() - 2] - 1);
		}
		else {
			weight = (GR_tr_sorted[pos + 1] - GR_tr_sorted[pos - 1]) / (GR_tr_count[pos + 1] + GR_tr_count[pos] + GR_tr_count[pos - 1] -1);
		}
		GR_tr_weights[i] = weight * 20;			
		if (weight<0){
			cout<<"Error: "<<weight<<endl;
			cout<<pos<<endl;
		}
	}
	//cout<<clone(GR_tr_weights).sort()<<endl;
 	cout<<"END"<<endl;
	return GR_tr_weights;
}


// [[Rcpp::depends(RcppEigen)]]
NumericMatrix MatrixMul(NumericMatrix x, NumericMatrix y) {
  Eigen::Map<Eigen::MatrixXd> X = as<Eigen::Map<Eigen::MatrixXd> >(x);
  Eigen::Map<Eigen::MatrixXd> Y = as<Eigen::Map<Eigen::MatrixXd> >(y);
  Eigen::MatrixXd Z = X * Y;
  return (Rcpp::NumericMatrix(wrap(Z)));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix MatrixMul_schur(NumericMatrix x, NumericMatrix y) { 
  arma::mat X = as<arma::mat>(x);
  arma::mat Y = as<arma::mat>(y); 
  arma::mat Z = X % Y;
  return(Rcpp::NumericMatrix(wrap(Z)));
} 
