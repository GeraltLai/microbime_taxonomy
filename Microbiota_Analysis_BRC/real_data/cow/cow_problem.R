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
tax = read.table(paste0(file,"cow_tax.txt"), header=TRUE)

### Metadata. Since we made this in Excel, not mothur, we can use the "row.names" modifier to automatically name the rows by the values in the first column (sample names)
meta = read.table(paste0(file,"cow_metadata.txt"), header=TRUE)

### SCFA data
SCFA = read.table(paste0(file,"cow_SCFA.txt"), header=TRUE)

### data
data=read.table(paste0(file,"cow_data.txt"),header=T)

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


#### Problem: Method 2 permutaion > sample size
new=data[data[,"family"]==names(table(data$family))[7],]
clean=select(OTU,row.names(new))
r=c()
for(i in 1:nrow(clean)){
  if(all(clean[i,]==0)) 
    r=c(r,i)
}
clean=clean[-r,];new_meta=meta[-r,]
dist=vegdist(clean, distance="bray")
adonis2(dist~AgeGroup,data=new_meta)

### last one row
new=data[data[,"family"]==names(table(data$family))[2],]
clean=select(OTU,row.names(new))
r=c()
for(i in 1:nrow(clean)){
  if(all(clean[i,]==0)) 
    r=c(r,i)
}
clean=clean[-r,];new_meta=meta[-r,]
dist=vegdist(clean, distance="bray")
adonis2(dist~AgeGroup,data=new_meta)

#### Problem: taxnomy
# different family but same genus
sum=c()
for(i in 1:length(table(data$family))){
  sum=c(sum,tax_f_g(names(table(data$family))[i]))
}
length(sum)
test=table(sum)
test1=subset(test,test>1)
test1
test2=subset(data,data$genus=="Clostridium")
table(test2$family)
