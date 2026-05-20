// PURPOSE:
//   Estimates the average probability of accepting a "bad" move
//   (a move that worsens the objective) at a given temperature T
//   in the simulated annealing (SA) framework.


//Ex. GetpBad(T=1)    → returns pBad=0.3  (too cold, try hotter)
//GetpBad(T=10)   → returns pBad=0.7  (still too cold)
//GetpBad(T=100)  → returns pBad=0.99 (too hot, go back)
//GetpBad(T=50)   → returns pBad=0.98 ✓ → this would be ideal for Tinit

//we want to be able to know what the probability of accepting a bad move is throughout all the temperatures. (not too high near the end; and not too low near the start)

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
double GetpBad(NumericMatrix G, NumericMatrix E, NumericMatrix X_tr,NumericMatrix Y_tr, NumericVector GR_tr, double SA_T, NumericVector GR_tr_weights, int obj_num);

//[[Rcpp::export]]
double GetpBad(NumericMatrix G, NumericMatrix E, NumericMatrix X_tr,NumericMatrix Y_tr, NumericVector GR_tr, double SA_T, NumericVector GR_tr_weights, int obj_num){

    //high objective is bad and means predictions are mostly wrong while obj == 1 means perfect predictions
	double obj = 1e30;
	int traits = G.nrow(), genes = G.ncol(), envi = E.ncol();
    	int SA_iter = 0, counter_bad = 0, N = 500000;
	int num_exp = GR_tr.length();
	double MAX_RATIO = 100;
	NumericVector circBuf;
	double sumCirc = 0;
//	int ran_tracer = 1;
	
	NumericMatrix TG_tr = MatrixMul(G, X_tr);
	NumericMatrix TE_tr = MatrixMul(E, Y_tr);

	float part = (float) genes / (genes + envi);

	NumericMatrix G_best = clone(G);
	NumericMatrix E_best = clone(E);
	NumericVector GR_guess_tr (GR_tr.length(), 0), temp_vector;
	char select, coerce_to_zero, deal_with_freedom;
	int whichRow, whichCol;

	double delta, sum_x, sum_y, sum_xy, sum_xx, sum_yy, pearson, GR_mean_diff, obj_new, temp_value, temp_value_2;
	double pBad = 1, SA_r, SA_s, SA_pAcc;
    	int freedomavail = num_exp, freedomusing = traits * (genes + envi);
 	float freedomratio, perturborzero;
//	cout<<"BEGIN"<<endl;
//	NumericVector GR_tr_sorted = clone(GR_tr).sort();
//	NumericVector GR_tr_count (unique(GR_tr_sorted).length(), 0);
//	//Rcout<<as<arma::rowvec>(GR_tr[seq_len(20)])<<endl;
//	//cout<<GR_tr_sorted<<endl;
//	NumericVector GR_tr_weights (GR_tr.length(), 0);
//	float weight;
//	int pos = 0;
//	double current = GR_tr_sorted[0];
//	
//	for (int i=0; i<GR_tr_sorted.length(); i++){
//		if (GR_tr_sorted[i] == current){
//			GR_tr_count[pos]++;
//		} else {
//			current = GR_tr_sorted[i];
//			pos++;
//			GR_tr_count[pos]++;
//		}	
//	}
//	GR_tr_sorted = unique(GR_tr_sorted).sort();
//	for (int i=0; i<GR_tr.length(); i++ ) {
//		pos = distance(GR_tr_sorted.begin(), std::find(GR_tr_sorted.begin(), GR_tr_sorted.end(), GR_tr[i]));
//		if (pos == 0){
//			weight = (GR_tr_sorted[1] - GR_tr_sorted[0]) / (GR_tr_count[0] + GR_tr_count[1] - 1);
//		}
//		else if (pos == GR_tr_sorted.length() - 1){
//			weight = (GR_tr_sorted[GR_tr_sorted.length() - 1] - GR_tr_sorted[GR_tr_sorted.length() - 2]) / (GR_tr_count[GR_tr_sorted.length() - 1] + GR_tr_count[GR_tr_sorted.length() - 2] - 1);
//		}
//		else {
//			weight = (GR_tr_sorted[pos + 1] - GR_tr_sorted[pos - 1]) / (GR_tr_count[pos + 1] + GR_tr_count[pos] + GR_tr_count[pos - 1] -1);
//		}
//		GR_tr_weights[i] = weight * 20;			
//		if (weight<0){
//			cout<<"Error: "<<weight<<endl;
//			cout<<pos<<endl;
//		}
//	}
//	//cout<<clone(GR_tr_weights).sort()<<endl;
// 	cout<<"END"<<endl;
	while (SA_iter < N){
		float current = (float) R::runif(0, 1);
		whichRow = int (floor(R::runif(0, traits)));
  		if (current <= part){
			select = 'G';
			whichCol = int (floor(R::runif(0, genes))); 
		} else {
			select = 'E';			
			whichCol = int (floor(R::runif(0, envi)));	
		}

		perturborzero = R::runif(0, 1);
		freedomratio = (float) freedomavail / (float) freedomusing;
		if (perturborzero <= freedomratio) {
			coerce_to_zero = 'N';
			delta = R::runif(-pBad, pBad);
			if (ISNAN(delta)){
				cout<<"Nan found in Delta!"<<endl;
				cout<<" Delta "<<delta<<" pBad "<<pBad<<" sumCirc "<<sumCirc<<endl;
				return 0;
			}
			if (select =='G'){
				if (G( whichRow, whichCol ) == 0){
					deal_with_freedom = '+';
				}
			        temp_value = G( whichRow, whichCol );			
				temp_vector = TG_tr( whichRow, _ );
				G( whichRow, whichCol ) += delta;
                if (G(whichRow, whichCol) > 1.0) {
                    G(whichRow, whichCol) = 1.0;
                }
                if (G(whichRow, whichCol) < -1.0) {
                    G(whichRow, whichCol) = -1.0;
                }
                delta = G(whichRow, whichCol) - temp_value;
     				TG_tr( whichRow, _ ) = TG_tr( whichRow, _ ) + delta * X_tr( whichCol, _ );
			} else {
				if (E( whichRow, whichCol ) == 0){
					deal_with_freedom = '+';
				}		
			        temp_value = E( whichRow, whichCol );			
				temp_vector = TE_tr( whichRow, _ );
				E( whichRow, whichCol ) += delta;
                if (E(whichRow, whichCol) > 1.0) {
                    E(whichRow, whichCol) = 1.0;
                }
                if (E(whichRow, whichCol) < -1.0) {
                    E(whichRow, whichCol) = -1.0;
                }
                delta = E(whichRow, whichCol) - temp_value;

                
     				TE_tr( whichRow, _ ) = TE_tr( whichRow, _ ) + delta * Y_tr( whichCol, _ );
    			}
		} else {
			coerce_to_zero = 'Y';
			if (select =='G'){
				if (G(whichRow, whichCol ) != 0){
					deal_with_freedom = '-';
					temp_value = G( whichRow, whichCol );
					G( whichRow, whichCol ) = 0;
					temp_vector = TG_tr( whichRow, _ );
     					TG_tr( whichRow, _ ) = TG_tr( whichRow, _ ) - temp_value * X_tr( whichCol, _ );
				} else {
					temp_value = 0;
					temp_vector = TG_tr( whichRow, _ );
				}
			} else {
				if (E( whichRow, whichCol ) != 0){
					deal_with_freedom = '-';
					temp_value = E ( whichRow, whichCol );
					E( whichRow, whichCol ) = 0;
					temp_vector = TE_tr( whichRow, _ );
     					TE_tr( whichRow, _ ) = TE_tr( whichRow, _ ) - temp_value * Y_tr( whichCol, _ );
    				} else {
					temp_value = 0;
					temp_vector = TE_tr( whichRow, _ );
				}	 
			}
		}
		NumericMatrix product = MatrixMul_schur(TG_tr, TE_tr);
		GR_mean_diff = 0;
 		for (int i=0; i<num_exp; i++) {
	         	GR_guess_tr[i] = sum(product( _, i )) / 1500;
			// Decides which objective function calculation to use 
            if (obj_num == 1) {
                temp_value_2 = max(GR_guess_tr[i], GR_tr[i]) / min(GR_guess_tr[i], GR_tr[i]);
                if (temp_value_2 > MAX_RATIO || temp_value_2 < 0) temp_value_2 = MAX_RATIO;
                GR_mean_diff += GR_tr_weights[i] * log(temp_value_2);
            } else if (obj_num == 2) {
                temp_value_2 = max(GR_guess_tr[i], GR_tr[i]) / min(GR_guess_tr[i], GR_tr[i]);
                if (temp_value_2 > MAX_RATIO || temp_value_2 < 0) temp_value_2 = MAX_RATIO;
                GR_mean_diff += GR_tr_weights[i] * temp_value_2;
            } else {
                temp_value_2 = abs(GR_guess_tr[i] - GR_tr[i]) / GR_tr[i];
                if (temp_value_2 > MAX_RATIO) temp_value_2 = MAX_RATIO;
                GR_mean_diff += GR_tr_weights[i] * temp_value_2;
            }
		}
		obj_new = (obj_num == 1) ? exp(GR_mean_diff / num_exp) : GR_mean_diff / num_exp; //decides based on objective
		SA_iter++;
		SA_pAcc = exp((obj-obj_new) / SA_T);
		SA_r = R::runif(0, 1);
		if (SA_r < SA_pAcc) {
			obj = obj_new;
			if (deal_with_freedom == '-'){
				freedomusing--;
			}
			if (deal_with_freedom == '+'){
				freedomusing++;
			}	
		} else {
			if (select == 'G'){
				G( whichRow, whichCol ) = temp_value;
     				TG_tr( whichRow, _ ) = temp_vector;
			} else {
				E( whichRow, whichCol ) = temp_value;
     				TE_tr( whichRow, _ ) = temp_vector;
    			}
			if (SA_pAcc < 0){
				cout<<"Warning: negative SA_pAcc!!"<<SA_pAcc<<'\n'<<endl;
				return 0;
			}
			if (sumCirc < 0){
				cout<<"Warning: negative sumCirc!!"<<sumCirc<<'\n'<<endl;
				return 0;
			}
			sumCirc += SA_pAcc;	
			counter_bad++;
		        pBad = sumCirc / (double) counter_bad;
		}
		deal_with_freedom = 'N';
	}
	cout<<pBad<<endl;
	if (counter_bad==0){
    		return(1);
  	}
  	else {
    		return(pBad);    
  	}
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
