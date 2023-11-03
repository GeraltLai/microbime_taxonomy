##### packages #####
library(ape)
library(dplyr)
library(ggplot2)
library(gplots)
library(lme4)
library(phangorn)
library(phyloseq)
library(plotly)
library(tidyr)
library(vegan)
library(VennDiagram)
library(tidyverse)
library(ecole) # zero-adjusted Bray–Curtis
library(gtools)  # permutations packages
library(ggupset) # upset plot
library(cowplot) # combine plot

###### read data
rm(list = ls())
data1 = readRDS(file = "C:/Users/lab205/Desktop/basic.rds")
data2 = readRDS(file = "C:/Users/lab205/Desktop/method3.rds")
data3 = readRDS(file = "C:/Users/lab205/Desktop/method4.rds")


### sort out different delta result
data = data3
total_genus = 6
#
delta_list = names(table(data[[1]]$delta1))
n = length(delta_list)
delta_name = sapply(1:n, function(n) paste0("delta=", delta_list[n]))
result_list = list()
total_case = rbind(data$genus1, data$genus2, data$genus3)
total_control = rbind(data$genus4, data$genus5, data$genus6)
for(i in 1:n){
  delta = delta_list[i]
  case = total_case[total_case$delta1 == delta, ]
  control = total_control[total_control$delta1 == delta, ]
  TP = case[case$individual == 1, ]
  FP = control[control$individual == 1, ]
  FN = case[case$individual == 0, ]
  TN = control[control$individual == 0, ]
  ### TP
  TP_with_non = table(factor(TP$with_nonsig, levels = 0:(total_genus-1)))
  TP_with_sig = table(factor(TP$with_sig, levels = 0:(total_genus-1)))
  TP_total = table(factor(TP$total, levels = 0:total_genus))
  new_TP = c(TP_with_non, TP_with_sig, TP_total, nrow(TP))
  ### FP
  FP_with_non = table(factor(FP$with_nonsig, levels = 0:(total_genus-1)))
  FP_with_sig = table(factor(FP$with_sig, levels = 0:(total_genus-1)))
  FP_total = table(factor(FP$total, levels = 0:total_genus))
  new_FP = c(FP_with_non, FP_with_sig, FP_total, nrow(FP))
  ### FN
  FN_with_non = table(factor(FN$with_nonsig, levels = 0:(total_genus-1)))
  FN_with_sig = table(factor(FN$with_sig, levels = 0:(total_genus-1)))
  FN_total = table(factor(FN$total, levels = 0:total_genus))
  new_FN = c(FN_with_non, FN_with_sig, FN_total, nrow(FN))
  ### TN
  TN_with_non = table(factor(TN$with_nonsig, levels = 0:(total_genus-1)))
  TN_with_sig = table(factor(TN$with_sig, levels = 0:(total_genus-1)))
  TN_total = table(factor(TN$total, levels = 0:total_genus))
  new_TN = c(TN_with_non, TN_with_sig, TN_total, nrow(TN))
  ### summary
  result = rbind(new_TP, new_FP, new_FN, new_TN)
  result_list[[i]] = result
}
names(result_list) = delta_name 
result_list[[1]]


### correct algorithm
data = result_list
total_genus = 6
#
output = list()
for(iii in 2:length(data)){
  list_name = names(data)[iii]
  correct_delta = data[[iii]]
  classification = c("with_nonsig", "with_sig", "total")
  with_nons_name = sapply(0:(total_genus-1), function(x) paste0("with_nonsig_", x))
  with_sig_name = sapply(0:(total_genus-1), function(x) paste0("with_sig_", x))
  total_count_name = sapply(0:total_genus, function(x) paste0("total_", x))
  colnames(correct_delta) = c(with_nons_name, with_sig_name, total_count_name, "sum")
  ###
  correct_algorithm_sig = list()
  list_id = 0
  classification_number = length(classification)
  for(ii in 1:classification_number){
    correct_delta_id = grep(classification[ii], colnames(correct_delta))
    correct_classification = correct_delta[,correct_delta_id]
    id = ncol(correct_classification)
    for(i in 1:(id-1)){
      list_id = list_id + 1
      keep_id = i
      classification_name = colnames(correct_classification)[i]
      TP_count = sum(correct_classification[1, 1:keep_id])
      FP_count = sum(correct_classification[2, 1:keep_id])
      FN_count = sum(correct_classification[3,]) + sum(correct_classification[1,-(1:keep_id)])
      TN_count = sum(correct_classification[4,]) + sum(correct_classification[2,-(1:keep_id)])
      correct_matrix = matrix(data = c(TP_count, FP_count, FN_count, TN_count),
                              ncol = 2, byrow = T)
      correct_matrix = as.data.frame(correct_matrix)
      colnames(correct_matrix) = c("truth+", "truth-")
      row.names(correct_matrix) = c("predict+", "predict-")
      correct_algorithm_sig[[list_id]] = correct_matrix
      names(correct_algorithm_sig)[[list_id]] = paste0("sig_result_keep_", classification_name)
    }
  }
  ###
  correct_algorithm_nonsig = list()
  list_id = 0
  classification_number = length(classification)
  for(ii in 1:classification_number){
    correct_delta_id = grep(classification[ii], colnames(correct_delta))
    correct_classification = correct_delta[,correct_delta_id]
    id = ncol(correct_classification)
    for(i in 1:(id-1)){
      list_id = list_id + 1
      keep_id = i
      classification_name = colnames(correct_classification)[i]
      TP_count = sum(correct_classification[1,]) + sum(correct_classification[3,-(1:keep_id)])
      FP_count = sum(correct_classification[2,]) + sum(correct_classification[4,-(1:keep_id)])
      FN_count = sum(correct_classification[3, 1:keep_id])
      TN_count = sum(correct_classification[4, 1:keep_id])
      correct_matrix = matrix(data = c(TP_count, FP_count, FN_count, TN_count),
                              ncol = 2, byrow = T)
      correct_matrix = as.data.frame(correct_matrix)
      colnames(correct_matrix) = c("truth+", "truth-")
      row.names(correct_matrix) = c("predict+", "predict-")
      correct_algorithm_nonsig[[list_id]] = correct_matrix
      names(correct_algorithm_nonsig)[[list_id]] = paste0("non_sig_result_keep_", classification_name)
    }
  }
  individual_output = c(correct_algorithm_sig, correct_algorithm_nonsig)
  output[[iii-1]] = individual_output
  names(output)[[iii-1]] = list_name
}

output[8]



