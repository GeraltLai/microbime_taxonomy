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
file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/cow/"

### OTU table
OTU = read.table(paste0(file,"cow_otu.txt"), header=TRUE)

### Taxonomy of each OTU
tax = read.table(paste0(file,"cow_tax.txt"),header=T)

### Metadata
meta = read.table(paste0(file,"cow_metadata.txt"),header=T)

### SCFA data
SCFA = read.table(paste0(file,"cow_SCFA.txt"),header=T)

##############################
##### data preprocessing #####
##############################
### separate tax data
col_name <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
temp_matrix <- matrix(nrow = NROW(tax), ncol = length(col_name),
                      dimnames = list(NULL, col_name))
sep1 <- strsplit(tax[["Taxonomy"]], split=";")
for(i in 1:NROW(tax)){
  each_taxa <- sep1[[i]]
  temp_matrix[i, ] <- substr(each_taxa, 4, nchar(each_taxa))
  print(i)
}
sep <- cbind(tax, temp_matrix)
rm(list=c("col_name","temp_matrix","sep1","each_taxa"))

### delete no OTU in sep data
data=c()
data <- lapply(names(OTU), function(x){
  index <- which(sep[,1] == x)
  data<- rbind(data, sep[index,])
})
data=do.call(rbind, data)
# write.table(data,file=paste0(file,"cow_data.txt"))
