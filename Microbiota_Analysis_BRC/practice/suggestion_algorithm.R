######
rm(list = ls())

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

#####################
##### load data #####
#####################
crop = "wildrice"
file = paste0("C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/",crop,"/")

### OTU table
OTU = read.table(paste0(file,crop,"_otu.txt"), header=T)

### Taxonomy of each OTU
tax = read.table(paste0(file,crop,"_tax.txt"), header=T)

### Metadata
meta = read.table(paste0(file,crop,"_metadata.txt"), header=T)

### data
data = read.table(paste0(file,crop,"_data.txt"), header=T)

### factor
label = "experimental_factor"

#### Group3
all_result = readRDS(file = paste0(file,"result/",crop,"_group_3",".rds"))

#### pairwise reult
pairwise_result = readRDS(file = paste0(file,"result/",crop,"_pairwise.rds"))

#### minusone result
minusone_result = readRDS(file = paste0(file,"result/",crop,"_minusone.rds"))

#####################
##### function ######
#####################
### pull genus result together by genus name
pull_result = function(all, pairwise, minus){
  genus = all[-1,]
  n = nrow(genus)
  row = nrow(genus)
  col = (row-1)*2+1
  table = matrix(data = 0, nrow = row, ncol = col)
  for(i in 1:n){
    step1 = genus[i,4]
    name = row.names(genus)[i]
    id = grep(name, row.names(pairwise))
    step2 = pairwise[id, 4]
    id2 = grep(name, row.names(minus))
    step3 = minus[id2, 4]
    output = c(step1, step2, step3)
    table[i,] = output
  }
  row.names(table)=row.names(genus)
  return(table)
}
### estimate_cutoff_value
estimate_cutoff_value = function(num_col, threshold) {
  value = threshold * num_col
  min_value = 1
  cutoff = max(value, min_value)
  return(cutoff)
}
### algorithm 1 : No weighted
algorithm_1 = function(data, threshold = 0) {
  num_col = ncol(data)
  result = matrix(0, nrow = nrow(data), ncol = ncol(data))
  index = data < 0.05
  result[index] = 1
  cutoff = estimate_cutoff_value(num_col, threshold)
  output = apply(result, 1, function(x) {
    sum(x) >= cutoff
  })
  output = as.data.frame(output)
  colnames(output) = "suggestion"
  row.names(output) = row.names(data)
  return(output)
}
### algorithm 2 : weighted
algorithm_2 = function(data, threshold = 0){
  n1 = nrow(data)
  result = matrix(0, nrow = nrow(data), ncol = ncol(data))
  index = data < 0.05
  result[index] = 1
  step1_weighted = n1
  step2_weighted = choose(n1, 2)
  step3_weighted = choose(n1,1)
  max_value = (1/step1_weighted) + 
    (n1 - 1)/step2_weighted+
    (n1 - 1)/step3_weighted
  threshold = 1
  value = threshold * max_value
  min_value = 1/step1_weighted
  value = max(c(value, min_value))
  output = c()
  for(i in 1:nrow(result)){
    step1 = 1
    step1_result = result[i, step1]/step1_weighted
    step2 = n1 -1 
    step2_region = step1 + 1:step2
    step2_result = sum(result[i, step2_region])/step2_weighted
    step3 = n1 -1 
    step3_region = (step1 + step2) + 1:step3
    step3_result = sum(result[i, step3_region])/step3_weighted
    part = step1_result + step2_result + step3_result
    part
    if(part >= value){
      output = rbind(output,"Yes")
    }else{
      output = rbind(output,"No")
    }
  }
  row.names(output) = row.names(data)
  colnames(output) = "suggestion"
  output = as.data.frame(output)
  return(output)
}

### test algorithm
i = 5
result1 = as.data.frame(all_result[[i]])
result2 = as.data.frame(pairwise_result[[i]])
result3 = as.data.frame(minusone_result[[i]])
test = pull_result(all = result1, pairwise = result2, minus = result3)
test
algorithm_1(data = test, threshold = 0.8)
algorithm_2(data = test, threshold = 0.9)

