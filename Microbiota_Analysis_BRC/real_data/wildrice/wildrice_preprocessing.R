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

##############################
##### data preprocessing #####
##############################
Size=lapply(1:ncol(OTU),function(x) sum(OTU[,x]))
Size=do.call(rbind,Size)
data=cbind(Size,tax)
data=data[,-ncol(data)]
data=data[complete.cases(data),]
# write.table(data,file=paste0(file,"wildrice_data.txt"))