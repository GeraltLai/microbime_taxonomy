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
# mu1 = 0.1
sigma1 = 1
# delta1 = 2
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


#saveRDS(OTU, file = paste0(file,"/otu.rds"))
#saveRDS(meta, file = paste0(file,"/metadata.rds"))
#saveRDS(tax, file = paste0(file,"/tax.rds"))



