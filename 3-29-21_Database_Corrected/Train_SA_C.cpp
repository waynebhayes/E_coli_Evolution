//OUTPUTS:
//   Returns a List containing:
//     [[1]] G          — optimized gene-trait matrix
//     [[2]] E          — optimized environment-trait matrix
//     [[3]] GR_guess_tr — predicted growth rates for training set
//     [[4]] SA_iter    — number of iterations taken to converge
//     [[5]] GR_mean_diff



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
List Train_SA_C(NumericMatrix G, NumericMatrix E, NumericMatrix X_tr,NumericMatrix Y_tr, NumericVector GR_tr, NumericMatrix X_va, NumericMatrix Y_va, NumericVector GR_va, int num_cycle, int num_round, int i, NumericVector T_range, NumericVector GR_tr_weights);

//[[Rcpp::export]]
List Train_SA_C(NumericMatrix G, NumericMatrix E, NumericMatrix X_tr,NumericMatrix Y_tr, NumericVector GR_tr, NumericMatrix X_va, NumericMatrix Y_va, NumericVector GR_va, int num_cycle, int num_round, int i, NumericVector T_range, NumericVector GR_tr_weights){
	cout<<"Firstone: "<<R::runif(0,1)<<endl;
	Environment global = Environment::global_env();
	double SA_numIters = (double) global["Iter"];
	double obj_best = global["obj_best"];
  	std::string filename_report = "run_report_" + to_string(i) + ".txt";
  	std::string filename_result = "run_result_" + to_string(i) + ".txt";
	cout<<"Cycle "<<num_cycle<<"Round "<<num_round<<endl;
 	
	ofstream outfile;
	outfile.open(filename_report, ios::out | ios::app);
	
	outfile<<"Cycle "<<num_cycle<<"Round "<<num_round<<endl;
   	outfile.close();
	int num_exp = GR_tr.length();
	double obj = 1e30, obj_prev = 1, obj_curr = 1;
	int traits = G.nrow(), genes = G.ncol(), envi = E.ncol();
    	int numCirc = 1000;
   	NumericVector circBuf;
	circBuf = rep(1, numCirc); 	
	double sumCirc = sum(circBuf);
	int iCirc = 0;
    	double pBad = 1;
	double MAX_RATIO = 100;
	cout<<"Reached T_range\n"<<endl;
	cout<<"Tinit is: "<<T_range[0]<<" Tfinal is: "<<T_range[1]<<endl;
    	double SA_Tinit = T_range[0];
   	double SA_Tfinal = T_range[1];

//	int ran_tracer = 1;
	int SA_iter = 0;
	double SA_s, SA_pAcc, SA_r;
	double SA_T = SA_Tinit;
	cout<<"Start computing SA_lambda"<<endl;
	double SA_lambda = -log(SA_Tfinal / SA_Tinit);
	cout<<"SA_lambda is "<<SA_lambda<<endl;
	int success = 0, failure = 0;
	NumericVector GR_guess_tr(GR_tr.length(), 0);
	NumericVector GR_guess_va(GR_va.length(), 0);
	cout<<"Finished establishing GR_guess_tr"<<endl;
	cout<<"Dim of G: "<<G.nrow()<<" x "<<G.ncol()<<endl;
	cout<<"Dim of E: "<<E.nrow()<<" x "<<E.ncol()<<endl;
	cout<<"Dim of X_tr: "<<X_tr.nrow()<<" x "<<X_tr.ncol()<<endl;
	cout<<"Dim of Y_tr: "<<Y_tr.nrow()<<" x "<<Y_tr.ncol()<<endl;
	NumericMatrix TG_tr = MatrixMul(G, X_tr);
	NumericMatrix TE_tr = MatrixMul(E, Y_tr);
	cout<<"Finished establishing TG_tr and TE_tr!"<<endl;
	float part = (float) genes / (genes + envi);
	NumericMatrix G_best = clone(G);
	NumericMatrix E_best = clone(E);
	NumericVector best = NumericVector::create(Named("iter",0), Named("obj")=0, Named("rho")=0, Named("GRdif")=0);
	NumericVector best_global, temp_vector;
	char select, coerce_to_zero, deal_with_freedom;
	int whichRow, whichCol;
	cout<<"Finished establishing best!"<<endl; 
	double delta, sum_x, sum_y, sum_xy, sum_xx, sum_yy, pearson, GR_mean_diff, obj_new, temp_value, temp_value_2;

	int status_freq = 500000;
    	int freedomavail = num_exp, freedomusing = traits * (genes + envi);
 	float freedomratio, perturborzero;
//	cout<<"BEGIN"<<endl;
//	NumericVector GR_tr_sorted = clone(GR_tr).sort();
//	NumericVector GR_tr_count (unique(GR_tr_sorted).length(), 0);
//	Rcout<<as<arma::rowvec>(GR_tr[seq_len(20)]`)<<endl;
//	cout<<GR_tr_sorted<<endl;
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
//	cout<<clone(GR_tr_weights).sort()<<endl;
// 	cout<<"END"<<endl;
 	cout<<"About to start!"<<endl;
	while (true){
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
			delta = R::runif(-pBad/10, pBad/10);
			if (select =='G'){
				if (G( whichRow, whichCol ) == 0){
					deal_with_freedom = '+';
				}
			        temp_value = G( whichRow, whichCol );			
				temp_vector = TG_tr( whichRow, _ );
				G( whichRow, whichCol ) += delta;
     				TG_tr( whichRow, _ ) = TG_tr( whichRow, _ ) + delta * X_tr( whichCol, _ );
			} else {
				if (E( whichRow, whichCol ) == 0){
					deal_with_freedom = '+';
				}		
			        temp_value = E( whichRow, whichCol );			
				temp_vector = TE_tr( whichRow, _ );
				E( whichRow, whichCol ) += delta;
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
//		if (SA_iter == 10000){
//			cout << "Address of temp_vector: " << &temp_vector << endl;
//			cout << "Address of TE_tr: " << &TE_tr << endl;
//			cout << "Address of TG_tr: " << &TG_tr << endl;
//			cout << "temp_vector[20] before:" << temp_vector[20] << endl;
//			cout << "TE_tr( whichRow, 20 ) before:" << TE_tr( whichRow, 20 ) << endl; 
//			cout << "TG_tr( whichRow, 20 ) before:" << TG_tr( whichRow, 20 ) << endl; 
//			cout << "Define temp_vector[20] = 200;" << endl;
//			temp_vector[20] = 200;
//			cout << "temp_vector[20] after:" << temp_vector[20] << endl;
//			cout << "TE_tr( whichRow, 20 ) after:" << TE_tr( whichRow, 20 ) << endl; 
//			cout << "TG_tr( whichRow, 20 ) after:" << TG_tr( whichRow, 20 ) << endl; 
//			cout << "Define TE_tr[whichRow, 20] = 800;" << endl;
//			TE_tr( whichRow, 20 ) = 800;
//			cout << "temp_vector[20] after:" << temp_vector[20] << endl;
//			cout << "TE_tr( whichRow, 20 ) after:" << TE_tr( whichRow, 20 ) << endl; 
//			cout << "TG_tr( whichRow, 20 ) after:" << TG_tr( whichRow, 20 ) << endl; 
//		}
	//	ran_tracer = ran_tracer+4;
	//	cout<<current<<"  "<<whichRow<<"  "<<whichCol<<"  "<<delta<<endl;
		//if (pBad > 1e-11) {
		//	delta = R::runif(-1, 1);
		//}
	 //	cout<<current<<'\t'<<part<<endl;
//		if (select =='G'){
//			G( whichRow, whichCol ) += delta;
//     			TG_tr( whichRow, _ ) = TG_tr( whichRow, _ ) + delta * X_tr( whichCol, _ );
//		} else {
//			E( whichRow, whichCol ) += delta;
 //    			TE_tr( whichRow, _ ) = TE_tr( whichRow, _ ) + delta * Y_tr( whichCol, _ );
//    		}
		NumericMatrix product = MatrixMul_schur(TG_tr, TE_tr);
		GR_mean_diff = 0;
 		for (int i=0; i<num_exp; i++) {
	         	GR_guess_tr[i] = sum(product( _, i )) / 1500;
			temp_value_2 = max(GR_guess_tr[i], GR_tr[i]) / min(GR_guess_tr[i], GR_tr[i]);
			if (temp_value_2 > MAX_RATIO || temp_value_2 < 0){
				temp_value_2 = MAX_RATIO;
			}
			GR_mean_diff += GR_tr_weights[i] * log(temp_value_2);
		}
		//GR_mean_diff = pow(GR_mean_diff, 1/num_exp);
//		sum_x = sum(GR_tr);
//		sum_y = sum(GR_guess_tr);
//		sum_xy = sum(GR_tr * GR_guess_tr);
//		sum_xx = sum(pow(GR_tr, 2));
//		sum_yy = sum(pow(GR_guess_tr, 2));
// 		pearson = (num_exp * sum_xy - sum_x * sum_y) / sqrt( (num_exp * sum_xx - sum_x * sum_x) * (num_exp * sum_yy - sum_y * sum_y) );
		
		
	//	obj_new = GR_mean_diff / (1 + pearson);
		obj_new = exp(GR_mean_diff / num_exp);
		SA_s = SA_iter / SA_numIters;
		SA_iter++;
		SA_T = SA_Tinit * exp (-SA_lambda * SA_s);
		SA_pAcc = exp((obj-obj_new) / SA_T);
		if (ISNAN(SA_pAcc)){
 			ofstream outfile;
			outfile.open(filename_report, ios::out | ios::app);
			outfile<<"Nan found in SA_pAcc!"<<endl;
			outfile<<"obj_new "<<obj_new<<" Delta "<<delta<<" pBad "<<pBad<<" pearson "<<pearson<<" sumCirc "<<sumCirc<<endl;
			outfile.close();
			return 0;
		}
		SA_r = R::runif(0, 1);
//		ran_tracer++;
		// cout<<SA_T<<endl;
		//if (SA_iter % 10000 == 0){
		//	cout<<ran_tracer<<'\t'<<endl;
	        //}	
		// cout<<SA_s<<'\t'<<SA_iter<<'\t'<<SA_numIters<<endl;
		//cout<<obj<<'\t'<<obj_new<<'\t'<<((obj-obj_new) / SA_T)<<'\t'<<SA_pAcc<<endl;
	        // cout<<SA_r<<"   "<<SA_pAcc<<obj_new<<SA_T<<endl;	
		if (SA_r < SA_pAcc) {
	  		if (obj > obj_new){
				G_best = clone(G);
        			E_best = clone(E);
       				best["iter"] = SA_iter;
        			best["obj"] = obj_new;
        			best["rho"] = pearson;
        			best["GRdif"] = GR_mean_diff;
				success++;
		      	}
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
			iCirc++;
			circBuf[iCirc % numCirc] = SA_pAcc;
			sumCirc = sum(circBuf);
			failure++;
		        pBad = sumCirc / (double) numCirc;
			if (SA_pAcc < 0){
				cout<<"Warning: negative SA_pAcc!!"<<SA_pAcc<<'\n'<<endl;
				return 0;
			}
			if (sumCirc < 0){
				cout<<"Warning: negative sumCirc!!"<<sumCirc<<'\n'<<endl;
				return 0;
			}
		}
		deal_with_freedom = 'N';
	//	if (abs(obj_curr - obj_prev) > 0.5){
	//		status_freq = 100;
	//	} else {
	//		status_freq = 10000;
	//	}
//		if (ran_tracer == 96666){
//			cout<<ran_tracer<<'\t'<<R::runif(0,1)<<endl;
//		}
		if (SA_iter % status_freq == 0) {
			cout<<"Iter "<<SA_iter<<" pBad "<<pBad<<" T "<<SA_T<<" Obj "<<obj<<" Rho "<<pearson<<" GRdif "<<GR_mean_diff<<" success "<<success<<" failure "<<failure<<" obj_new "<<obj_new<<" SA_pAcc "<<SA_pAcc<<" Delta "<<delta<<" SA_r "<<SA_r<<" freedomratio "<<freedomratio<<" i: "<<i<<endl;
//			cout<<"tracer: "<<ran_tracer<<'\t'<<"random value: "<<R::runif(0,1)<<endl;
 			ofstream outfile;
			outfile.open(filename_report, ios::out | ios::app);
			outfile<<"Iter "<<SA_iter<<" pBad "<<pBad<<" T "<<SA_T<<" Obj "<<obj<<" Rho "<<pearson<<" GRdif "<<GR_mean_diff<<" success "<<success<<" failure "<<failure<<" obj_new "<<obj_new<<" SA_pAcc "<<SA_pAcc<<" Delta "<<delta<<" SA_r "<<SA_r<<" freedomratio "<<freedomratio<<endl;
   			outfile.close();
			if (SA_iter >= SA_numIters) {
				cout<<"Cannot move further"<<endl;
	 			cout<<best<<endl;
				NumericMatrix pred_current = MatrixMul_schur(MatrixMul(G_best, X_tr), MatrixMul(E_best, Y_tr)) / 1500;
 				NumericVector pred_rate(pred_current.ncol());
				for (int i = 0; i < pred_current.ncol(); i++) {
					pred_rate[i] = sum(pred_current( _, i ));
				}
				NumericMatrix pred_current_va = MatrixMul_schur(MatrixMul(G_best, X_va), MatrixMul(E_best, Y_va)) / 1500;
 				NumericVector pred_rate_va(pred_current_va.ncol());
				for (int i = 0; i < pred_current_va.ncol(); i++) {
					pred_rate_va[i] = sum(pred_current_va( _, i ));
				}
				if (best[1] < best_global[0]) {
					best_global = best[Rcpp::Range(1,3)];
					global["best_global"] = best_global;
					cout<<"best vector"<<best_global<<endl;
				}
				ofstream outfile;
				outfile.open(filename_report, ios::out | ios::app);
				outfile<<"Cannot move further"<<endl;
				outfile<<"iter "<<best[0]<<" obj"<<best[1]<<" rho "<<best[2]<<" GRdif "<<best[3]<<endl;
				outfile<<"Training:"<<'\n'<<"Lab Growth Rates:"<<'\n'<<GR_tr<<endl;
				outfile<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<pred_rate<<endl;
				outfile<<"Testing:"<<'\n'<<"Lab Growth Rates:"<<'\n'<<GR_va<<endl;
				outfile<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<pred_rate_va<<endl;
				outfile<<"Best performance so far: "<<'\n'<<"Obj: "<<best_global[0]<<"Rho: "<<best_global[1]<<"Diff: "<<best_global[2]<<endl;
   				outfile.close();
				
				cout<<best<<endl;
				cout<<"Train: "<<'\n'<<"Lab Growth Rates:"<<'\n'<<GR_tr<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<pred_rate<<endl;	
				cout<<"Testing: "<<'\n'<<"Lab Growth Rates:"<<'\n'<<GR_va<<endl;
				cout<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<pred_rate_va<<endl;
				num_cycle++;
				std::scientific;	
				ofstream outfile_2;
		 		outfile_2.open(filename_result, ios::out | ios::app);
				outfile_2<<"G after training: "<<G_best<<'\n'<<"E after training: "<<E_best<<endl;
   				outfile_2.close();
				std::cout.precision(7);
				List runagain;
				if (num_cycle>3) {
					num_cycle = 1;
					num_round++;

					if (num_round > 3) {
						cout<<"Cannot get a good solution. Raise the goal of obj by 100%"<<endl;
						obj_best = obj_best * 2;
						num_round = 1;
						cout<<"Now the goal of obj is: "<<obj_best<<'\n'<<"Round "<<num_round<<"Cycle "<<num_cycle<<endl;

						ofstream outfile;
						outfile.open(filename_report, ios::out | ios::app);
						outfile<<"Cannot get a good solution. Raise the goal of obj by 100%"<<'\n';
						outfile<<"Now the goal of obj is: "<<obj_best<<'\n'<<"Round "<<num_round<<"Cycle "<<num_cycle<<endl;
   						outfile.close();
					}

						cout<<"G, E redefined. Anneal again with new G and E as input"<<endl;
						
						ofstream outfile;
		 				outfile.open(filename_report, ios::out | ios::app);
						outfile<<"G, E redefined. Anneal again with new G and E as input"<<endl;
   						outfile.close();
					
						NumericVector G_initial = rnorm(traits * genes);
						G_initial.attr("dim") = Dimension(traits, genes);
						global["G_initial"] = G_initial;
						
						NumericVector E_initial = rnorm(traits * envi);
						E_initial.attr("dim") = Dimension(traits, envi);
						global["E_initial"] = E_initial;

						G = as<NumericMatrix>(G_initial);	
						E = as<NumericMatrix>(E_initial);
			
						ofstream outfile_2;
		 				outfile_2.open(filename_result, ios::out | ios::app);
						outfile_2<<"Initial G: "<<G<<'\n'<<"Initial E: "<<E<<endl;
   						outfile_2.close();

 						runagain = Train_SA_C(G, E, X_tr, Y_tr, GR_tr, X_va, Y_va, GR_va, num_cycle, num_round, i, T_range, GR_tr_weights);
				} else{
					cout<<"Anneal again with the best G and E so far as input"<<endl;							
					ofstream outfile;
		 			outfile.open(filename_report, ios::out | ios::app);
					outfile<<"Anneal again with the best G and E so far as input"<<endl;
   					outfile.close();

					runagain = Train_SA_C(G, E, X_tr, Y_tr, GR_tr, X_va, Y_va, GR_va, num_cycle, num_round, i, T_range, GR_tr_weights);
				}
				return(runagain);	
			}
		}
//		if (iter == 5870000) {
//			status_freq = 100;
//		}	
		if (obj < obj_best | (num_cycle==2 & SA_iter == SA_numIters-1)) {

			NumericMatrix pred_current_va = MatrixMul_schur(MatrixMul(G_best, X_va), MatrixMul(E_best, Y_va)) / 1500;
			for (int i = 0; i < pred_current_va.ncol(); i++) {
				GR_guess_va[i] = sum(pred_current_va( _, i ));
			}

			cout<<"Good enough"<<'\n';
			cout<<"iter "<<SA_iter<<"obj "<<obj<<"rho "<<pearson<<"GRdif "<<GR_mean_diff<<'\n';
			cout<<"Train: "<<'\n'<<"Lab Growth Rates:"<<'\n'<<GR_tr<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<GR_guess_tr<<endl;	
		
			ofstream outfile;
			outfile.open(filename_report, ios::out | ios::app);
			outfile<<"Good enough"<<'\n';
			outfile<<"iter "<<best[0]<<" obj"<<best[1]<<" rho "<<best[2]<<" GRdif "<<best[3]<<endl;
			outfile<<"Training:"<<'\n'<<"Lab Growth Rates:"<<'\n'<<GR_tr<<endl;
			outfile<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<GR_guess_tr<<endl;
			outfile<<"Testing:"<<'\n'<<"Lab Growth Rates:"<<'\n'<<GR_va<<endl;
			outfile<<"Dot Production Prediction of Lab Growth Rates:"<<'\n'<<GR_guess_va<<endl;
			outfile<<"Best performance so far: "<<'\n'<<"Obj: "<<best_global[0]<<"Rho: "<<best_global[1]<<"Diff: "<<best_global[2]<<endl;
   			outfile.close();

			ofstream outfile_2;
		 	outfile_2.open(filename_result, ios::out | ios::app);
			outfile_2<<"G after training: "<<G_best<<'\n'<<"E after training: "<<E_best<<endl;
   			outfile_2.close();
			return (List::create(G, E, GR_guess_tr, SA_iter, GR_mean_diff));		
		}				
	}
}

// [[Rcpp::depends(RcppEigen)]]
NumericMatrix MatrixMul(NumericMatrix x, NumericMatrix y) {
  Eigen::Map<Eigen::MatrixXd> X = as<Eigen::Map<Eigen::MatrixXd> >(x);
  Eigen::Map<Eigen::MatrixXd> Y = as<Eigen::Map<Eigen::MatrixXd> >(y);
  Eigen::MatrixXd Z = X * Y;
  return (wrap(Z));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix MatrixMul_schur(NumericMatrix x, NumericMatrix y) { 
  arma::mat X = as<arma::mat>(x);
  arma::mat Y = as<arma::mat>(y); 
  arma::mat Z = X % Y;
  return(Rcpp::NumericMatrix(wrap(Z)));
} 
