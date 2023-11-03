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

#########################
##### set parameter #####
#########################
### data
file = paste0("C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/simulation_data")
group = 2
sample = 4
total_sample = sample * group
genus= 6
genus_OTU = 3
total_OTU = genus * genus_OTU
summary_result = list()

### parameter
mu1 = 1
sigma1 = 1
delta1 = 2
delta0 = 0
case_mu2 = mu1 + delta1
control_mu2 = mu1 + delta0
sigma2 = sigma1

### log normal simulation data
data = rlnorm(sample * total_OTU, meanlog = mu1, sdlog = sigma1)
data_matrix = matrix(data = data, nrow = sample , ncol = total_OTU, byrow = T)
case_data = rlnorm((sample * total_OTU)/2, meanlog = case_mu2, sdlog = sigma2)
case_matrix = matrix(data = case_data, nrow = sample, ncol = total_OTU/2 )
control_data = rlnorm((sample * total_OTU)/2, meanlog = control_mu2, sdlog = sigma2)
control_matrix = matrix(data = control_data, nrow = sample, ncol = total_OTU/2 )
OTU = cbind(case_matrix, control_matrix)
OTU = rbind(data_matrix, OTU)
row_name = c()
for(i in 1:nrow(OTU)){
  row_name = c(row_name, paste0('sample',i))
}
OTU_name = c()
for(i in 1:ncol(OTU)){
  OTU_name = c(OTU_name, paste0('OTU',i))
}
colnames(OTU) = OTU_name
row.names(OTU) = row_name
OTU = as.data.frame(round(OTU,0))
index = OTU == 0
OTU[index] = 1

### meta data simulation
factor = rep(letters[1:group], each = sample)
meta = matrix(data = factor, nrow = group * sample)
colnames(meta) = "factor"
row.names(meta) = row_name
meta = as.data.frame(meta)

### tax data simulation
tax_family_name = rep("family", total_OTU)
tax_genus_name = c()
for(i in 1:genus){
  genus_name = rep(paste0('genus', i), genus_OTU)
  tax_genus_name = c(tax_genus_name,  genus_name)
}
tax = data.frame(tax_family_name,tax_genus_name )
row.names(tax) = OTU_name
colnames(tax) = c("family", "genus")



##################################################
###### step1 : all genus individual output #######
##################################################
step1_output = c()
dist = bray0(OTU)
adonis2(dist ~ factor, data = meta)
result = adonis2(dist ~ factor, data = meta)[1,3:5]
size = round(sum(OTU),0)
result = cbind(size, result)
step1_output = rbind(step1_output, result)
genus = unique(tax$genus)
for(i in genus){
  genus_group = subset(tax, tax[["genus"]] == i)
  genus_id = row.names(genus_group)
  new_OTU = select(OTU, all_of(genus_id))
  size = round(sum(new_OTU), 0)
  dist = bray0(new_OTU)
  result = adonis2(dist ~ factor, data = meta)[1,3:5]
  result = cbind(size, result)
  step1_output = rbind(step1_output, result)
}
row.names(step1_output) = c("total", genus)
colnames(step1_output) = c("size", "R2", "F", "Pr(>F)")
step1_output = as.data.frame(step1_output)
step1_output


############################################
###### step2 : pairwise genus output #######
############################################
limit = 0
pick = step1_output
pick1 = subset(pick, pick[["size"]] >= max(pick[["size"]]) * limit)
pick1 = pick1[-1,]
genus = unique(tax$genus)
step2_output = c()
name = c()
for(i in 1:(nrow(pick1)-1))
  for(j in (i+1):nrow(pick1)){
    size = sum(pick1[i,1] + pick1[j,1])
    genus_id1 = genus[i]
    genus_id2 = genus[j]
    name = c(name, paste0(genus_id1,"-",genus_id2))
    genus_group1 = subset(tax, tax[["genus"]] == genus_id1)
    genus_group2 = subset(tax, tax[["genus"]] == genus_id2)
    genus_OTU = c(row.names(genus_group1),row.names(genus_group2))
    clean = select(OTU, all_of(genus_OTU))
    dist = bray0(clean)
    ans = adonis2(dist ~ factor, data = meta)[1,3:5]
    output = round(c(size, ans[1,1], ans[1,2], ans[1,3]), 4)
    step2_output = rbind(step2_output, output)
  }
colnames(step2_output) = c("size", "R2", "F", "Pr(>F)")
row.names(step2_output) = name
step2_output = as.data.frame(step2_output)
step2_output

#############################################
###### step3 : minus one genus output #######
#############################################
limit = 0
step3_output = c()
name = c()
pick = step1_output
pick1 = subset(pick, pick[["size"]] >= max(pick[["size"]]) * limit)
pick1 = pick1[-1,]
all_gene = row.names(pick1)
for(i in row.names(pick1)){
  delet = which(row.names(pick1) == i)
  size = sum(pick1[-delet, 1])
  part_gene = all_gene[-(delet)]
  name = c(name, paste0(part_gene, collapse = "-"))
  part_tax = filter(tax, tax[["genus"]]!= i)
  part_OTU = row.names(part_tax)
  clean = select(OTU, all_of(part_OTU))
  dist = bray0(clean)
  ans = adonis2(dist ~ factor, data = meta)[1,3:5]
  output = round(c(size, ans[1,1], ans[1,2], ans[1,3]), 4)
  step3_output = rbind(step3_output, output)
}
colnames(step3_output) = c("size","R2","F","Pr(>F)")
row.names(step3_output) = name
step3_output = as.data.frame(step3_output)
step3_output

###########################################
###### suggest algorithm function #########
###########################################
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
### count table
count_table = function(data){
  result = matrix(0, nrow = nrow(data), ncol = ncol(data))
  index = data < 0.05
  result[index] = 1
  output = apply(result, 1, function(x) {
    sum(x)
  })
  output = as.data.frame(output)
  colnames(output) = "count"
  row.names(output) = row.names(data)
  return(output)
}
### algorithm 2 : weighted table
algorithm_2 = function(data){
  n1 = nrow(data)
  result = matrix(0, nrow = nrow(data), ncol = ncol(data))
  index = data < 0.05
  result[index] = 1
  step1_weighted = n1
  step2_weighted = choose(n1, 2)
  step3_weighted = choose(n1,1)
  output = c()
  for(i in 1:nrow(result)){
    step1 = 1
    step1_result = round(result[i, step1]/step1_weighted,3)
    step2 = n1 -1 
    step2_region = step1 + 1:step2
    step2_result = round(sum(result[i, step2_region])/step2_weighted,3)
    step3 = n1 -1 
    step3_region = (step1 + step2) + 1:step3
    step3_result = round(sum(result[i, step3_region])/step3_weighted,3)
    part = step1_result + step2_result + step3_result
    output = rbind(output,part)
  }
  row.names(output) = row.names(data)
  colnames(output) = "weighted_count"
  output = as.data.frame(output)
  return(output)
}

### compare result
test1 = pull_result(all = step1_output, pairwise = step2_output, minus = step3_output)
individual_result = step1_output[,4]
my_result = count_table(test1)
weighted_result = algorithm_2(test1)
if(individual_result[1] >= 0.05){
  my_result = rbind(0, my_result)
  weighted_result = rbind(0, weighted_result)
}else{
  my_result = rbind(1, my_result)
  weighted_result = rbind(1, weighted_result)
}
row.names(my_result)[1] = "family"
row.names(weighted_result)[1] = "family"
compare_result = cbind(individual_result, my_result, weighted_result)
list_name = paste0("mu=",mu1,"|","delta=",delta1)
summary_result[[list_name]] = compare_result
step1_output
summary_result



