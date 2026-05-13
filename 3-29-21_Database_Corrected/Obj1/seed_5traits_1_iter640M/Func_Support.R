Shuffle<-function(exp_index, paperdiv){
  
  exp_index_shuffle<-exp_index
  num_papers<-length(unique(exp_index[ ,2]))
  
  if (paperdiv){
    
    for (i in 1:num_papers){
      exp_paper<-matrix(exp_index[exp_index[ ,2]==i, ], ncol=2)
      exp_paper_shuffle<-exp_paper[sample(nrow(exp_paper)), ]
      exp_index_shuffle[exp_index_shuffle[ ,2]==i, ]<-exp_paper_shuffle
    }
  }
  
  else {
    exp_index_shuffle<-exp_index_shuffle[sample(exp_index_shuffle[ ,1]), ]
  }
  
  return (exp_index_shuffle)
}

Kfold<-function(seed, K, paperdiv, sample_ID){
  
  set.seed(seed)
  
  num_exp<-length(sample_ID)
  
  dash_pos<-regexpr("-", sample_ID)
  attributes(dash_pos)<-NULL
  
  paper_ID<-rep(0, num_exp)
  
  for (i in 1:num_exp){
    paper_ID[i]<-substr(sample_ID[i], 1, dash_pos[i]-1)
  }
  
  papers<-unique(paper_ID)
  
  num_papers<-length(papers)
  
  from_paper<-rep(0, num_exp)
  
  for (i in 1:num_exp){
    from_paper[i]<-grep(paper_ID[i], papers)
  }
  
  exp_index<-cbind(1:num_exp, from_paper)
  
  exp_index_shuffle<-Shuffle(exp_index, paperdiv)
  
  num_papers<-length(unique(exp_index_shuffle[ ,2]))
  num_validation_samples<-vector("integer")
  
  tr_row_index<-list()
  va_row_index<-list()
  
  for (i in 1:K){
    tr_row_index[[i]]<-vector("integer")
    va_row_index[[i]]<-vector("integer")
  }
  
  if (paperdiv){
    rows_acc<-0
    for (i in 1:num_papers){
      all_exp<-matrix(exp_index[exp_index[,2]==i, ], ncol=2)
      num_selected<-floor(nrow(all_exp)/K)
      num_validation_samples<-append(num_validation_samples, num_selected)
      
      for (j in 1:K){
        if (num_selected!=0){
          row_index<-(num_selected*(j-1)+1):(num_selected*j)
          va_row_index[[j]]<-append(va_row_index[[j]], row_index+rows_acc)
        }
        else{
          row_index<-NULL
        }
        
        tr_row_index[[j]]<-append(tr_row_index[[j]], setdiff(1:nrow(all_exp), row_index)+rows_acc)
      }
      rows_acc<-rows_acc+nrow(all_exp)
    }
  }  
  else {
    for (j in 1:K){
      num_validation_samples<-floor(length(sample_ID/K))
      row_index<-(num_validation_samples*(j-1)+1):(num_validation_samples*j)
      va_row_index[[j]]<-append(va_row_index[[j]], row_index)
      tr_row_index[[j]]<-append(tr_row_index[[j]], setdiff(1:length(sample_ID), row_index))
    }
  }
  
  return (list(tr_row_index, va_row_index, exp_index_shuffle[ ,1]))
}

Validate<-function(G, E, X_va, Y_va, GR_va){
  TG_va<-G%*%X_va
  TE_va<-E%*%Y_va
  
  GR_guess_va<-rep(0, length(GR_va))
  
  for (i in 1:length(GR_guess_va)){
    # GR_guess[i]<-(TG[,i]+L[[3]])%*%TE[,i]
    GR_guess_va[i]<-TG_va[,i]%*%TE_va[,i]/1500
  }
  
  rho_sp <- cor(GR_va, GR_guess_va, method = "spearman")
  rho_pr <- cor(GR_va, GR_guess_va, method = "pearson")
  
  GR_mean_diff<-mean(abs(GR_guess_va-GR_va))
  
  rho <- (0*rho_sp + 2*rho_pr)/2
  
  obj<-abs(GR_mean_diff*(1-rho))
  
  return(list(obj, GR_guess_va))
  
}

Draw<-function(K, result_input, keyword){
  
  progress<-c("0%", paste(round(100*(1:K)/K), "%", sep=""))
  
  filename<-paste(keyword, ".jpeg", sep="")
  
  jpeg(file=filename, width=240*2*K, height=240*2*K)
  layout(matrix(c(1:K, 1:K),K,2))
  step<-length(result_input[[1]])/K
  j<-0
  
  for (i in 1:K){
    
    lab_data<-result_input$Lab[((i-1)*step+1):(i*step)]
    simu_data<-result_input$Simu[((i-1)*step+1):(i*step)]
    
    order_data<-order(lab_data)
    lab_data_ordered<-lab_data[order_data]
    simu_data_rearr<-simu_data[order_data]
    
    plot(lab_data_ordered, main=paste(progress[i], progress[i+1], sep="-"), ylim=c(min(c(lab_data, simu_data)), max(c(lab_data, simu_data))), cex.main=K/2, cex.lab=K/2, cex.axis=K/2, cex=K/2)
    points(simu_data_rearr, col=2, cex=K/2)
    legend("bottomright", legend=c("Real Data", "Simulation"), pch=c(1,1), col=c("black", "red"))
  }
  
  dev.off()
  
  return()
}

None_zero_idx_same<-function(m1, m2){
  
  
  product_multi<-m2*m1
  product_minus<-m2-m1
  
  result<-which(product_multi!=0, arr.ind=T)
  
  value_change<-rep(0, nrow(result))
  sign_change<-rep("?", nrow(result))
  
  for (i in 1:nrow(result)){
    idx<-as.vector(result[i, ])
    r<-idx[1]
    c<-idx[2]
    value_change[i]<-product_minus[r,c]
    sign_change[i]=case_when(
      (sign(m1[r,c])==-1 && sign(m2[r,c])==-1) ~ "-",
      (sign(m1[r,c])==1 && sign(m2[r,c])==1) ~ "+",
      (sign(m1[r,c])==-1 && sign(m2[r,c])==1) ~ "- -> +",
      (sign(m1[r,c])==1 && sign(m2[r,c])==-1) ~ "+ -> -"
    )
    if (is.na(sign_change[i])) {
      print(m1[idx])
      print(m2[idx])
    }
  }
  
  
  result<-cbind(result, sign_change, value_change)
  
  colnames(result)<-c("Row", "Col", "Sign Change", "Value Change")
  
  return(result)
}

Delete<-function(DB, Colname){
  
  Num<-grep(Colname, colnames(DB))
  newtable_secondhalf<-DB[,(Num+1):ncol(DB)]
  newtable_firsthalf<-DB[,1:(Num-1)]
  
  name_secondhalf<-colnames(DB)[(Num+1):ncol(DB)]
  name_firsthalf<-colnames(DB)[1:(Num-1)]
  
  newDB<-cbind(newtable_firsthalf, newtable_secondhalf)
  colnames(newDB)<-c(name_firsthalf, name_secondhalf)
  
  return(newDB)
}

Ins_right<-function(DB, Num, NewName){
  
  New<-rep("", nrow(DB))
  
  if (Num < ncol(DB)){
    table_secondhalf<-DB[,(Num+1):ncol(DB)]
    name_firsthalf<-colnames(DB)[1:Num]
    name_secondhalf<-colnames(DB)[(Num+1):ncol(DB)]
    
    DB<-cbind(DB, New)
    
    DB[, Num+1]<-New
    DB[, (Num+2):ncol(DB)]<-table_secondhalf
    colnames(DB)<-c(name_firsthalf, NewName, name_secondhalf)
  }
  
  else {
    name_new<-c(colnames(DB), NewName)
    DB<-cbind(DB, New)
    colnames(DB)<-name_new
  }
  return(DB)
}

Preprocessing<-function(input, synthesis = F){
  
  if (is.character(input)){
    DB<-read.csv(input, stringsAsFactors = FALSE)[,-1]  
  } else {
    DB<-input
  }
  
  
  colnames_total<-colnames(DB)
  col_gene<-(grep("Gene_start", colnames_total)+1):(grep("Growth_Rate", colnames_total)[1]-1)
  col_gene_start<-grep("Gene_start", colnames_total)+1
  
  if (synthesis==T){
    DB_wild<-data.frame()
    
    for (i in 1:nrow(DB)){
      if (sum(DB[i, col_gene])==0){
        DB_wild<-rbind(DB_wild, DB[i,])
      }
    }
    
    return(DB_wild)
  }
  
  i<-col_gene_start
  empty_col<-rep(0, nrow(DB))
  while(i<grep("Growth_Rate", colnames(DB))[1]){
    
    if (grepl("_Deletion", colnames(DB)[i])){
      underscore_pos<-regexpr("_Deletion", colnames(DB)[i])
      attributes(underscore_pos)<-NULL
      gene_name<-substr(colnames(DB)[i], 1, underscore_pos-1)
      gene_all_con<-grep(gene_name, colnames(DB))
      
      new_col<-empty_col
      
      for (j in 1:nrow(DB)){
        if (sum(DB[j, gene_all_con])==0){
          new_col[j]<-1
        }
      }
      
      DB<-Ins_right(DB, i, paste(gene_name, "_Normal", sep=""))
      DB[ ,paste(gene_name, "_Normal", sep="")]<-new_col
      
      if (sum(DB[, i])==0){
        DB<-Delete(DB, colnames(DB)[i])
        i<-i+1
      }
      else{
        i<-i+2
      }
      next;
    }
    
    if (sum(DB[, i])==0){
      DB<-Delete(DB, colnames(DB)[i])
    }
    else{
      i<-i+1
    }
  }
  
  if (is.character(input)){
    write.csv(DB, "Data_Matrix_Reduced.csv")  
  } else {
    write.csv(DB, "Data_Matrix_Reduced_wild.csv")
  }
  
  return(DB)
}
