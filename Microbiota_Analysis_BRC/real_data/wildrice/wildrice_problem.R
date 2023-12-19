#####
rm(list=ls())

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
#####################
##### load data #####
#####################
file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/wildrice/"
### OTU table
OTU = read.table(paste0(file,"wildrice_otu.txt"),header=T)

### Taxonomy of each OTU
tax = read.table(paste0(file,"wildrice_tax.txt"),header=T)

### Metadata
meta = read.table(paste0(file,"wildrice_metadata.txt"), header=TRUE)

### data
data=read.table(paste0(file,"wildrice_data.txt"),header=T)

#######################
###### function #######
#######################
### tax family to genus function
tax_f_g=function(x){
  new=subset(data,data$family==x)
  return(names(table(new$genus)))
}
### tax genus to family function
tax_g_f=function(x){
  new=subset(data,data$genus==x)
  return(names(table(new$family)))
}


#### Problem: taxanomy
sum=c()
for(i in 1:length(table(data$family))){
  sum=c(sum,tax_f_g(names(table(data$family))[i]))
}
length(sum)
test=table(sum)
test1=subset(test,test>1)
test1
test2=subset(data,data$genus=="Incertae Sedis")
table(test2$family)
