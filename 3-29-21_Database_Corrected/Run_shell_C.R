library(stringi)
library(stringr)
library(dplyr)
library(RcppEigen)
library(RcppArmadillo)
library(Rcpp)
source("../Main_new_C_par.R")
source("../Func_Support.R")
args=commandArgs(T)  #creates a character vector captured by the command line
input<-args[1]
seed<-as.integer(args[2])
traits<-as.integer(args[3])
obj_best<<-as.numeric(args[4])
Iter<<-as.integer(args[5])
remaining<-as.integer(args[6:length(args)])
#Iter<<-Iter*traits

result<-Main_new_C_par(input, seed, traits, obj_best, Iter, remaining)
