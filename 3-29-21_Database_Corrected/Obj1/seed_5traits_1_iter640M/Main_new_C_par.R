Main_new_C_par<-function(input, seed, traits, obj_best, Iter, exp_info){

  library(foreach)
  library(doParallel)
  library(doMC)
  DB<-read.csv(input, fileEncoding="latin1")
  set.seed(seed)    
  val_result<-list("Lab"=c(), "Simu"=c())
  tr_result<-list("Lab"=c(), "Simu"=c())
  G_each_K<-c()
  E_each_K<-c()
  va_idx_all<-c()
  tr_idx_all<-c()
  
  # colnames_total<-colnames(DB)
  # col_gene<-(grep("Gene_start", colnames_total)+1):(grep("Growth_Rate", colnames_total)-1)
  # col_gene_start<-grep("Gene_start", colnames_total)+1
  # col_env<-(grep("Env_start", colnames_total)+1):(grep("Gene_start", colnames_total)[1]-1)
  # 
    

  # DB_env<<-DB[ , col_env]
  # DB_gene<<-DB[ , col_gene]
  GR<-DB[, "Growth.rate..1.h."]
  # GR_basal<-DB[, "Basal_Growth_Rate"]
  # X<-t(DB_gene)
  # X[is.na(X)]<-0
  # Y<-t(DB_env)
  # Y[is.na(Y)]<-0
  
  Strain_info<-DB[, (match("Strain_Start",colnames(DB))+1):(match("Strain_End",colnames(DB))-1)]
  Medium_info<-DB[, (match("Medium_Start",colnames(DB))+1):(match("Medium_End",colnames(DB))-1)]
  Treatment_info<-DB[, (match("Treatment_Start",colnames(DB))+1):(match("Treatment_End",colnames(DB))-1)]
  
  X<-DB[,(match("Gene_Start",colnames(DB))+1):(match("Gene_End",colnames(DB))-1)]
  Y<-cbind(Strain_info, Medium_info, Treatment_info)
  X<-t(X)
  Y<-t(Y)
  X[is.na(X)]<-0
  Y[is.na(Y)]<-0
  print(dim(X))
  print(dim(Y))
  # if (basal){
  #   GR<-GR_lab-GR_basal  
  # }
  # else{
  #   GR<-GR_lab
  # }
  
  # paperdiv<-TRUE
  # 
  # output<-Kfold(seed, K, paperdiv, sample_ID)
  # tr_whichrow<-output[[1]]
  # va_whichrow<-output[[2]]
  # exp_index_shuffle<-output[[3]]
  
  num_papers<-length(unique(as.character(DB[, "Author"])))
  papers<-unique(DB[, "Author"])
  #print(papers) 
  #print(num_papers)
  tr_row_index<-list()
  va_row_index<-list()
  
  K<-length(exp_info)/3
  for (i in 1:K){
    print(exp_info[((i-1)*3+2)]:exp_info[((i-1)*3+3)])
    for_va<-exp_info[((i-1)*3+2)]:exp_info[((i-1)*3+3)]
    if (length(for_va)==1){
      va_row_index[[i]]<-for_va 
    } else{
      va_row_index[[i]]<-sample(for_va)
    }
    tr_row_index[[i]]<-sample((1:nrow(DB))[-for_va])
    print(for_va)	
  }
  G_initial<-list()
  E_initial<-list()
  
  
  result<-list()
  
   
  score<-vector("numeric")
  validation_score<-vector("numeric")
  
  GR_va_guess<-c()
  GR_va_all<-c()
  result_train<- vector(mode = "list", length = K)
  ls(.GlobalEnv)
  G_initial<-list()
  E_initial<-list()
  X_tr<-list()
  Y_tr<-list()
  X_va<-list()
  Y_va<-list()
  GR_tr<-list()
  GR_va<-list()
  weight<-list()
  
  
  cat("Sourcing Getweight.cpp\n")    
  sourceCpp("./Getweight.cpp") 
  cat("Sourcing Finished\n")
  for (i in 1:K){
    current_exp<-exp_info[(i-1)*3+1] 
    G_initial[[i]]<-matrix(rnorm(traits*nrow(X)), traits, nrow(X))
    E_initial[[i]]<-matrix(rnorm(traits*nrow(Y)), traits, nrow(Y))
    G<-G_initial[[i]]
    E<-E_initial[[i]]
    
    filename_report<-paste("run_report_", current_exp, ".txt", sep="")
    filename_result<-paste("run_result_", current_exp, ".txt", sep="")
    cat("seed = ", seed, "traits = ", traits, "\n", file=filename_report)
    cat("seed = ", seed, "traits = ", traits, "\n", file=filename_result)
   
    cat("Training: part ", current_exp, "\n")
    cat("Training: part ", current_exp, "\n", file=filename_report, append=T)
    
    
    cat("Initial G:\n", G_initial[[i]], "\n", file=filename_result, append=T)
    cat("Initial E:\n", E_initial[[i]], "\n", file=filename_result, append=T)
    
    
    cat("Training: part ", current_exp, "\n", file=filename_report, append=T)
    

    tr_index<-tr_row_index[[i]]   #<-------variable outside of loop
    va_index<-va_row_index[[i]]   #<-------variable outside of loop
    
    
    X_tr[[i]]<-X[, tr_row_index[[i]]]           #<-------variable outside of loop
    Y_tr[[i]]<-Y[, tr_row_index[[i]]]           #<-------variable outside of loop
    X_va[[i]]<-X[, va_row_index[[i]]]
    Y_va[[i]]<-Y[, va_row_index[[i]]]
    if (length(va_row_index[[i]])==1){
      X_va[[i]]<-matrix(X_va[[i]], ncol=1)
      Y_va[[i]]<-matrix(Y_va[[i]], ncol=1)
    }

    GR_tr[[i]]<-GR[tr_row_index[[i]]]           #<-------variable outside of loop
    GR_va[[i]]<-GR[va_row_index[[i]]]

    cat(length(GR_tr[[i]]), "experiments used for training \n")
    cat(length(GR_va[[i]]), "experiments used for testing, they are ", va_row_index[[i]], "\n") 
    cat("Round ", 1, "\n", "Goal of obj ", obj_best, "\n", file=filename_report, append=T)
    cat("Round ", 1, "\n", "Goal of obj ", obj_best, "\n")
    	
    weight[[i]] <- Getweight(GR_tr[[i]])
 }
  cat("Sourcing GetpBad.cpp\n")    
  sourceCpp("./GetpBad.cpp") 
  cat("Sourcing Finished\n")
  print(weight[[1]])  
  threshold<-0.98
  SA_Tinit<-1
  registerDoMC(detectCores(logical=F))
  TT<-foreach (i = 1:K) %dopar%{
    current_exp<-exp_info[(i-1)*3+1] 
    filename_report<-paste("run_report_", current_exp, ".txt", sep="")
    filename_result<-paste("run_result_", current_exp, ".txt", sep="")
    registerDoMC(detectCores(logical=F))
    while(GetpBad(G_initial[[i]], E_initial[[i]], X_tr[[i]], Y_tr[[i]], GR_tr[[i]], SA_Tinit, weight[[i]])<threshold){
       	SA_Tinit<-SA_Tinit*10
    }
    SA_Tinit<-SA_Tinit/2
    while(GetpBad(G_initial[[i]], E_initial[[i]], X_tr[[i]], Y_tr[[i]], GR_tr[[i]], SA_Tinit, weight[[i]])>threshold){
    	SA_Tinit<-SA_Tinit/2
    }
    cat("Experiment:", current_exp, "has Tinit of: ", SA_Tinit, '\n')     
    SA_Tfinal<-SA_Tinit/10
    
    threshold<-1e-6
    while(GetpBad(G_initial[[i]], E_initial[[i]], X_tr[[i]], Y_tr[[i]], GR_tr[[i]], SA_Tfinal, weight[[i]])>threshold){      
    	SA_Tfinal<-SA_Tfinal/10
    }
    SA_Tfinal<-SA_Tfinal*2
    while(GetpBad(G_initial[[i]], E_initial[[i]], X_tr[[i]], Y_tr[[i]], GR_tr[[i]], SA_Tfinal, weight[[i]])<threshold){
    	SA_Tfinal<-SA_Tfinal*2
    }
    cat("Experiment:", current_exp, "has Tfinal of: ", SA_Tfinal, '\n')     
#   T_range[[i]]<-c(500, 1e-5)
    T_range<-c(SA_Tinit, SA_Tfinal)
    cat("Temperature range: ", T_range, "\n")  
    cat("Temperature range: ", T_range, "\n", file=filename_report, append=T)
    T_range
  }  
  print("Rcpp compilation started")
  cat("Rcpp version: ", packageDescription("Rcpp")$Version, "\n")
  sourceCpp("./Train_SA_C.cpp") 
  print("Compilation finished")
  cat("Iter is: ", Iter, "\n")
    
#  foreach(i=1:K, .export=c("Validate", "GetpBad", "sourceCpp", "traits", "X", "Y", "tr_row_index", "va_row_index", "result_train", "G_each_K", "E_each_K", "tr_result", "Iter", "obj_best", "seed", "K")) %dopar% {
  re_<-foreach(i=1:K) %dopar% {

    	current_exp<-exp_info[(i-1)*3+1] 
        result<-Train_SA_C(G_initial[[i]], E_initial[[i]], X_tr[[i]], Y_tr[[i]], GR_tr[[i]], X_va[[i]], Y_va[[i]], GR_va[[i]], 1, 1, current_exp, TT[[i]], weight[[i]])
	print(length(result))
	print(class(result))
	return(result)
  }
  stopImplicitCluster()
  print(length(re_))  
  print(length(re_[[1]]))
    #score<-append(score, result_train[[i]][[1]])
  for (i in 1:K){ 
    current_exp<-exp_info[(i-1)*3+1] 
    filename_report<-paste("run_report_", current_exp, ".txt", sep="")
    filename_result<-paste("run_result_", current_exp, ".txt", sep="")
    G<-re_[[i]][[1]]
    cat("G=", t(G), "\n", file=filename_result, append=T)
    
    
    #G_prime<-G
    #G_prime[abs(G_prime)<sort(abs(G_prime), decreasing=TRUE)[length(sample_ID)]]<-0
    #cat("G_prime:\n", t(G_prime), "\n", file="run_result_top.txt", append=T)
    
    #if (i>1){
    #  G_prime_same_idx<-None_zero_idx_same(G_prime_prev, G_prime)
    #  cat("Non-zero element positions where G_prime shares with the previous one:\n", file="run_result_top.txt", append=T)
    #  write.table(G_prime_same_idx,  file="run_result_top.txt", append=T)
    
    #}
    #G_prime_prev<-G_prime
    
    
    
    E<-re_[[i]][[2]]
    cat("E=", t(E), "\n", file=filename_result, append=T)
    #E_prime<-E
    # if (length(sample_ID)<nrow(E_prime)*ncol(E_prime)){
    #  E_prime[abs(E_prime)<sort(abs(E_prime), decreasing=TRUE)[length(sample_ID)]]<-0
    # }
    
    
    # cat("E_prime:\n", E_prime, "\n", file="run_result_top.txt", append=T)
    # if (i>1){
    #   E_prime_same_idx<-None_zero_idx_same(E_prime_prev, E_prime)
    #   cat("Non-zero element positions where E_prime shares with the previous one:\n", file="run_result_top.txt", append=T)
    #   write.table(E_prime_same_idx,  file="run_result_top.txt", append=T)
    #   
    # }
    # E_prime_prev<-E_prime
    # 
    GR_tr_guess<-re_[[i]][[3]]
    
    tr_result$Lab<-append(tr_result$Lab, GR_tr[[i]])
    tr_result$Simu<-append(tr_result$Simu, GR_tr_guess)
    
    G_each_K<-append(G_each_K, G)
    E_each_K<-append(E_each_K, E)

    result_validate<-Validate(G, E, X_va[[i]], Y_va[[i]], GR_va[[i]])
    
    #validation_score<-append(validation_score, result_validate[[1]])
    
    val_result$Lab<-append(val_result$Lab, GR_va[[i]])
    val_result$Simu<-append(val_result$Simu, result_validate[[2]])
    
    cat("Test:\n", "Lab Growth Rates:\n")
    cat("Test:\n", "Lab Growth Rates:\n", file=filename_report, append=T)
    print(GR_va[[i]])
    cat(GR_va[[i]], file=filename_report, append=T)
    cat("Dot Product Prediction of Growth Rate:\n")
    cat("Dot Product Prediction of Growth Rate:\n", file=filename_report, append=T)
    cat(result_validate[[2]], file=filename_report, append=T)
  }  
  
  Draw(K, tr_result, "training")
  Draw(K, val_result, "testing")
  print("Finished")
  
  return(list(score, validation_score, tr_result, val_result, G_each_K, E_each_K, tr_row_index, va_row_index))
}
