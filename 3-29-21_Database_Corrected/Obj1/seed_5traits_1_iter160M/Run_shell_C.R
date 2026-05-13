my_paths<-.libPaths() 
my_paths<-c(my_paths, "/home/weir4/R/x86_64-pc-linux-gnu-library/3.5", "/pkg/R/3.5.1-9.25.2018/lib64/R/library", "/home/weir4/local/R_libs")
.libPaths(my_paths)
#
#library(stringi, lib.loc="~/local/R_libs")
#library(stringr, lib.loc="~/local/R_libs")
library(stringr)
library(dplyr)
library(RcppEigen)
library(RcppArmadillo)
library(Rcpp)
source("Main_new_C_par.R")
source("Func_Support.R")
args=commandArgs(T)
input<-args[1]
seed<-as.integer(args[2])
traits<-as.integer(args[3])
obj_best<<-as.numeric(args[4])
Iter<<-as.numeric(args[5])
remaining<-as.integer(args[6:length(args)])


result<-Main_new_C_par(input, seed, traits, obj_best, Iter, remaining)
