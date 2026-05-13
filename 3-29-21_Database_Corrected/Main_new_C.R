Main_new_C<-function(input, seed, traits, K, obj_best, Iter){

  DB<-read.csv(input, fileEncoding="latin1")
  DB_ori<-DB
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
  
  Strain_info<-DB[, match("Strain_Start",colnames(DB)):match("Strain_End",colnames(DB))]
  Medium_info<-DB[, match("Medium_Start",colnames(DB)):match("Medium_End",colnames(DB))]
  Treatment_info<-DB[, match("Treatment_Start",colnames(DB)):match("Treatment_End",colnames(DB))]
  
  X<-DB[,match("Gene_Start",colnames(DB)):match("Gene_End",colnames(DB))]
  Y<-cbind(Strain_info, Medium_info, Treatment_info)
  X<-t(X)
  Y<-t(Y)
  X[is.na(X)]<-0
  Y[is.na(Y)]<-0
  exp_keep<-rep(TRUE, nrow(DB))
  condi_keep_X<-rep(TRUE, nrow(X))
  for (i in 1:nrow(X)){
    non_zero<-which((X[i,]!=0))
    if (length(non_zero) <3){
      exp_keep[non_zero]<-FALSE
      condi_keep_X[i]<-FALSE      
    }
  }
  
  condi_keep_Y<-rep(TRUE, nrow(Y))
  for (i in 1:nrow(Y)){
    non_zero<-which(Y[i,]!=0)
    if (length(non_zero) <3){
      exp_keep[non_zero]<-FALSE
      condi_keep_Y[i]<-FALSE      
   }
  }
  # print(head(rownames(Y)[condi_keep_X]))
  X<-X[condi_keep_X, exp_keep]
  Y<-Y[condi_keep_Y, exp_keep]
  X<-rbind(X, rep(1, ncol(X)))
  Y<-rbind(Y, rep(1, ncol(Y)))
  DB<-DB[exp_keep, ]
  print(dim(X))
  print(dim(Y))
  write.csv(DB, "Merged_Data_4.csv")
  write.csv(X, "X.csv")
  write.csv(Y, "Y.csv")
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
  
  for (i in 1:K){
    tr_row_index[[i]]<-vector("integer")
    va_row_index[[i]]<-vector("integer")
  }
  for(i in 1:num_papers){
    # print(papers[i])
    exp_idx<-c()
    paper_search<-DB[, "Author"]==papers[i]
    for (k in 1:length(paper_search)){
      if (paper_search[k]){
        exp_idx<-append(exp_idx, k)
      }
    }
    num_selected<-floor(length(exp_idx)/K)
   # print(sort(exp_idx))
    for (j in 1:K){
      if (num_selected==0){
        tr_row_index[[j]]<-append(tr_row_index[[j]], exp_idx)
        next;
      }
       
      if (j!=K){
        row_selected<-exp_idx[(num_selected*(j-1)+1):(num_selected*j)]
        va_row_index[[j]]<-append(va_row_index[[j]], row_selected)
      }
      
      else {
        row_selected<-exp_idx[(num_selected*(j-1)+1):length(exp_idx)]
        va_row_index[[j]]<-append(va_row_index[[j]], row_selected)
      }
         
      tr_row_index[[j]]<-append(tr_row_index[[j]], setdiff(exp_idx, row_selected))

    
    }
  
   # print(length(tr_row_index[[1]]))
  }
  
  for (i in 1:K){
    
    tr_row_index[[i]]<-sample(tr_row_index[[i]])
    va_row_index[[i]]<-sample(va_row_index[[i]])
  }
  
  result<-list()
  
   
  score<-vector("numeric")
  validation_score<-vector("numeric")
  
  GR_va_guess<-c()
  GR_va_all<-c()
  result_train<- vector(mode = "list", length = K)
  cat("seed = ", seed, "traits = ", traits, "K = ", K, "\n", file="run_report.txt")
  cat("seed = ", seed, "traits = ", traits, "K = ", K, "\n", file="run_result.txt")
  
  for (i in 1:K){
    
    G_initial<<-matrix(rnorm(traits*nrow(X)), traits, nrow(X))
    E_initial<<-matrix(rnorm(traits*nrow(Y)), traits, nrow(Y))
    
    G<-G_initial
    E<-E_initial
   
    cat("Training: part ", i, " in ", K, "\n")
    cat("Training: part ", i, " in ", K, "\n", file="run_report.txt", append=T)
    
    
    cat("Initial G:\n", G, "\n", file="run_result.txt", append=T)
    cat("Initial E:\n", E, "\n", file="run_result.txt", append=T)
    
    
    cat("Training: part ", i, " in ", K, "\n", file="run_result.txt", append=T)
    cat("Training: part ", i, " in ", K, "\n", file="run_result_top.txt", append=T)
    
    best_global <<- c(10,10,10)
    names(best_global) <<- c("Obj", "Rho", "Diff") 
    tr_index<-tr_row_index[[i]]
    va_index<-va_row_index[[i]]
    
    
    X_tr<-X[, tr_index]
    Y_tr<-Y[, tr_index]
    X_va<-X[, va_index]
    Y_va<-Y[, va_index]
    GR_tr<-GR[tr_index]
    GR_va<-GR[va_index]
    cat(length(GR_tr), "experiments used for training", "\n")
    cat(length(GR_va), "experiments used for testing", "\n") 
    SA_Tinit<-1
    
    num_cycle<<-1
    num_round<<-1
    
    cat("Round ", num_round, "\n", "Goal of obj ", obj_best, "\n", file="run_report.txt", append=T)
    cat("Round ", num_round, "\n", "Goal of obj ", obj_best, "\n")
    
#    threshold<-0.5
#    while(GetpBad(G, E, X_tr, Y_tr, GR_tr, SA_Tinit, threshold)<threshold){
#      SA_Tinit<-SA_Tinit*10
#    }
#    SA_Tinit<-SA_Tinit/2
#    while(GetpBad(G, E, X_tr, Y_tr, GR_tr, SA_Tinit, threshold)>threshold){
#     SA_Tinit<-SA_Tinit/2
#    }
#    
#    SA_Tfinal<-SA_Tinit/10
#   
#    threshold<-1e-10
#    while(GetpBad(G, E, X_tr, Y_tr, GR_tr, SA_Tfinal, threshold)>threshold){      
#     SA_Tfinal<-SA_Tfinal/10
#    }
#    SA_Tfinal<-SA_Tfinal*2
#    while(GetpBad(G, E, X_tr, Y_tr, GR_tr, SA_Tfinal, threshold)<threshold){
#     SA_Tfinal<-SA_Tfinal*2
#    }
    T_range<-c(500, 1e-5)
#    T_range<-c(SA_Tinit, SA_Tfinal)
    cat("Temperature range: ", T_range, "\n")  
    cat("Temperature range: ", T_range, "\n", file="run_report.txt", append=T)
    
    print("Rcpp compilation started")
    cat("Rcpp version: ", packageDescription("Rcpp")$Version, "\n")
    sourceCpp("./Train_SA_C.cpp") 
    print("Compilation finished")
    cat("Iter is: ", Iter, "\n")
    result_train[[i]]<-Train_SA_C(G, E, X_tr, Y_tr, GR_tr, X_va, Y_va, GR_va, obj_best, T_range, Iter)
    
    score<-append(score, result_train[[i]][[1]])
    
    G<-result_train[[i]][[1]]
    cat("G=", t(G), "\n", file="run_result.txt", append=T)
    
    
    #G_prime<-G
    #G_prime[abs(G_prime)<sort(abs(G_prime), decreasing=TRUE)[length(sample_ID)]]<-0
    #cat("G_prime:\n", t(G_prime), "\n", file="run_result_top.txt", append=T)
    
    #if (i>1){
    #  G_prime_same_idx<-None_zero_idx_same(G_prime_prev, G_prime)
    #  cat("Non-zero element positions where G_prime shares with the previous one:\n", file="run_result_top.txt", append=T)
    #  write.table(G_prime_same_idx,  file="run_result_top.txt", append=T)
    
    #}
    #G_prime_prev<-G_prime
    
    
    
    E<-result_train[[i]][[2]]
    cat("E=", t(E), "\n", file="run_result.txt", append=T)
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
    GR_tr_guess<-result_train[[i]][[5]]
    
    tr_result$Lab<-append(tr_result$Lab, GR_tr)
    tr_result$Simu<-append(tr_result$Simu, GR_tr_guess)
    
    G_each_K<-append(G_each_K, G)
    E_each_K<-append(E_each_K, E)
    print(class(G))
    print(class(X_va))
    
    result_validate<-Validate(G, E, X_va, Y_va, GR_va)
    
    validation_score<-append(validation_score, result_validate[[1]])
    
    val_result$Lab<-append(val_result$Lab, GR_va)
    val_result$Simu<-append(val_result$Simu, result_validate[[2]])
    
    cat("Test:\n", "Lab Growth Rates:\n")
    cat("Test:\n", "Lab Growth Rates:\n", file="run_report.txt", append=T)
    print(GR_va)
    cat(GR_va, file="run_report.txt", append=T)
    cat("Dot Product Prediction of Growth Rate:\n")
    cat("Dot Product Prediction of Growth Rate:\n", file="run_report.txt", append=T)
    print(result_validate[[2]])
    cat(result_validate[[2]], file="run_report.txt", append=T)
    
  }
  Draw(K, tr_result, "training")
  Draw(K, val_result, "testing")
  
  return(list(score, validation_score, tr_result, val_result, G_each_K, E_each_K, tr_row_index, va_row_index))
}
