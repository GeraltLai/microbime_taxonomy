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


#####
step2_output = c()
name = c()
length = length(unique(tax$genus))
for(i in 1:(length-1))
  for(j in (i+1):length){
    genus_id1 = unique(tax$genus)[i]
    genus_id2 = unique(tax$genus)[j]
    name = c(name, paste0(genus_id1,"-",genus_id2))
    genus_group1 = subset(tax, tax[["genus"]] == genus_id1)
    genus_group2 = subset(tax, tax[["genus"]] == genus_id2)
    genus_OTU = c(row.names(genus_group1),row.names(genus_group2))
    clean = select(OTU, all_of(genus_OTU))
    size = sum(clean)
    ####################################
    dist = diag(x = 1, nrow = nrow(clean))
    for(ii in 1:(nrow(clean)-1))
      for(jj in (ii+1):nrow(clean)){
        A = clean[ii,] ; B = clean[jj,]
        top = sapply(1:length(A), function(x) (A[x]/sum(A)) * (B[x]/sum(B)))
        down1 = sapply(1:length(A), function(x) (A[x]/(sum(A))) ^ 2)
        down2 = sapply(1:length(A), function(x) (B[x]/(sum(B))) ^ 2)
        top = do.call(cbind, top) ; down1 = do.call(cbind, down1) ; down2 = do.call(cbind, down2)
        dist[ii,jj] = sum(top)/(sqrt(sum(down1) * sum(down2)))
      }
    label = t(dist)[lower.tri(dist)]
    # After transposing, remove the triangle and fill the lower triangle
    dist[lower.tri(dist)]=t(dist)[lower.tri(dist)]
    dist = as.dist(dist)
    dist = 1 - dist
    #####################################
    ans = adonis2(dist ~ factor, data = meta)[1,3:5]
    output = round(c(size, ans[1,1], ans[1,2], ans[1,3]), 4)
    step2_output = rbind(step2_output, output)
  }
colnames(step2_output) = c("size", "R2", "F", "Pr(>F)")
row.names(step2_output) = name
step2_output = as.data.frame(step2_output)
step2_output

##########
step2_output = c()
name = c()
######
i = unique(tax$genus)[1]
j = unique(tax$genus)[2]
genus_id1 = i
genus_id2 = j
name = c(name, paste0(genus_id1,"-",genus_id2))
genus_group1 = subset(tax, tax[["genus"]] == genus_id1)
genus_group2 = subset(tax, tax[["genus"]] == genus_id2)
genus_OTU = c(row.names(genus_group1),row.names(genus_group2))
clean = select(OTU, all_of(genus_OTU))
size = sum(clean)
####################################
dist = diag(x = 1, nrow = nrow(clean))
for(ii in 1:(nrow(clean)-1))
  for(jj in (ii+1):nrow(clean)){
    A = clean[ii,] ; B = clean[jj,]
    top = sapply(1:length(A), function(x) (A[x]/sum(A)) * (B[x]/sum(B)))
    down1 = sapply(1:length(A), function(x) (A[x]/(sum(A))) ^ 2)
    down2 = sapply(1:length(A), function(x) (B[x]/(sum(B))) ^ 2)
    top = do.call(cbind, top) ; down1 = do.call(cbind, down1) ; down2 = do.call(cbind, down2)
    dist[ii,jj] = sum(top)/(sqrt(sum(down1) * sum(down2)))
  }
label = t(dist)[lower.tri(dist)]
# After transposing, remove the triangle and fill the lower triangle
dist[lower.tri(dist)]=t(dist)[lower.tri(dist)]
dist = as.dist(dist)
dist = 1 - dist
#####################################
ans = adonis2(dist ~ factor, data = meta)
output = round(c(size, ans[1,1], ans[1,2], ans[1,3]), 4)
step2_output = rbind(step2_output, output)
##################
colnames(step2_output) = c("size", "R2", "F", "Pr(>F)")
row.names(step2_output) = name
step2_output = as.data.frame(step2_output)
step2_output



##### SIMI
A = c(1,10,14,10)
B = c(100,95,16,20)
top = sapply(1:length(A), function(x) (A[x]/sum(A)) * (B[x]/sum(B)))
down1 = sapply(1:length(A), function(x) (A[x]/(sum(A))) ^ 2)
down2 = sapply(1:length(A), function(x) (B[x]/(sum(B))) ^ 2)
sum(top)/(sqrt(sum(down1) * sum(down2)))




##### SIMI test
A = c(350,250,575,300,25)
B = c(350,275,275,300,300)
top = sapply(1:length(A), function(x) (A[x]/sum(A)) * (B[x]/sum(B)))
down1 = sapply(1:length(A), function(x) (A[x]/(sum(A))) ^ 2)
down2 = sapply(1:length(A), function(x) (B[x]/(sum(B))) ^ 2)
sum(top)/(sqrt(sum(down1) * sum(down2)))






