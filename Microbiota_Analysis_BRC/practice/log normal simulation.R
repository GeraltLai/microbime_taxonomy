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
set.seed(110)
group = 2
sample = 4
total_OTU = 4
mu1 = 1
sigma1 = 1
mu2 = mu1 * 3
sigma2 = sigma1 * 2

### log normal simulation data
data1 = rlnorm(sample * total_OTU, meanlog = mu1, sdlog = sigma1)
data2 = rlnorm(sample * total_OTU, meanlog = mu2, sdlog = sigma2)
OTU = matrix(data = c(data1, data2), 
             ncol = total_OTU, nrow = sample * group, byrow = T)
row_name = c()
for(i in 1:nrow(OTU)){
  row_name = c(row_name, paste0('sample',i))
}
col_name = c()
for(i in 1:ncol(OTU)){
  col_name = c(col_name, paste0('OTU',i))
}
colnames(OTU) = col_name
row.names(OTU) = row_name
OTU = as.data.frame(OTU)
OTU

### meta data simulation
factor = rep(letters[1:group], each = sample)
meta = matrix(data = factor, nrow = group * sample)
colnames(meta) = "factor"
row.names(meta) = row_name
meta = as.data.frame(meta)
meta

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
for(i in 1:ncol(OTU)){
  new_OTU = OTU[,i]
  size = round(sum(new_OTU), 0)
  dist = bray0(new_OTU)
  result = adonis2(dist ~ factor, data = meta)[1,3:5]
  result = cbind(size, result)
  step1_output = rbind(step1_output, result)
}
row.names(step1_output) = c("total", colnames(OTU))
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
step2_output = c()
name = c()
for(i in 1:(nrow(pick1)-1))
  for(j in (i+1):nrow(pick1)){
    size = sum(pick1[i,1] + pick1[j,1])
    id1 = rownames(pick1)[i]
    id2 = rownames(pick1)[j]
    name = c(name, paste0(id1,"-",id2))
    clean = OTU[, colnames(OTU) == id1 | colnames(OTU) == id2]
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
  clean = OTU[, colnames(OTU) != i]
  dist = bray0(clean)
  ans = adonis2(dist ~ factor, data = meta)[1,3:5]
  output = round(c(size, ans[1,1], ans[1,2], ans[1,3]), 4)
  step3_output = rbind(step3_output, output)
}
colnames(step3_output) = c("size","R2","F","Pr(>F)")
row.names(step3_output) = name
step3_output = as.data.frame(step3_output)
step3_output

###########################
###### ggupset plot #######
###########################
total = rbind(step1_output, step2_output, step3_output)
total = total[-1,]
total
genus = step1_output[-1,]
genus_name = row.names(genus)
ggplot(total, aes(x = reorder(row.names(total), `Pr(>F)`), y= `Pr(>F)`)) +
  geom_bar(stat = "identity") +
  axis_combmatrix(sep = "-", levels = c(genus_name)) +
  geom_hline(aes(yintercept=0.05),linetype='dashed', colour = "blue")+
  xlab("Genus combination")+
  ylab("p-value")

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

###
test1 = pull_result(all = step1_output, pairwise = step2_output, minus = step3_output)
test1
algorithm_1(data = test1, threshold = 0.8)
algorithm_2(data = test1, threshold = 0.9)
