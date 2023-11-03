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
crop="wildrice"
file=paste0("C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/",crop,"/")

### OTU table
OTU = read.table(paste0(file,crop,"_otu.txt"), header=T)

### Taxonomy of each OTU
tax = read.table(paste0(file,crop,"_tax.txt"),header=T)

### Metadata
meta = read.table(paste0(file,crop,"_metadata.txt"),header=T)

### data
data = read.table(paste0(file,crop,"_data.txt"),header=T)

### factor
label="experimental_factor"

#######################
###### function #######
#######################
### Bray–Curtis distance matrix 
dm=function(x, method,factor){
  if(length(unique(ts))==length(unique(data$family))){
    new=data[data[["family"]]==unique(x),]
  }else if(length(unique(ts))==length(unique(data$genus))){
    new=data[data[["genus"]]==unique(x),]
  }
  clean=select(OTU,row.names(new))
  if(method==1){
    dist=as.matrix(vegdist(clean, distance="bray"))
    index=is.nan(dist)
    dist[index]=1
    return(adonis2(dist~meta[, factor])[1,3:5])
  }
  if(method==2){
    index=lapply(1:nrow(clean),function(x){
      if(all(clean[x,]==0))
        index=c(x)
    })
    index=do.call(cbind, index)
    clean=clean ; new_meta=meta
    if(length(index)>1){
      clean=clean[-index,];new_meta=meta[-index,]
    }
    if((length(table(new_meta[[factor]]))<=1)){
      return(matrix(rep(0,3),nrow=1,ncol=3))
    }else{
      dist=vegdist(clean, distance="bray")
      return(adonis2(dist~new_meta[, factor])[1,3:5])
    }
  }
  if(method==3){
    dist=bray0(clean)
    return(adonis2(dist~meta[, factor])[1,3:5])
  }
}
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

##############################
##### family result ##########
##############################
family=list()
ts=data$family
for(a in 1:3){
  all_result=list()
  print(a)
  for(ff in unique(ts)){
    print(ff)
    all_result[[ff]]=dm(ff,a,label)
  }
  family[a]=list(all_result)
}
names(family)=c("method1","method2","method3")
rm(list=c("all_result"))

#############################
##### genus result ##########
#############################
genus=list()
ts=data$genus
for(a in 1:3){
  all_result=list()
  for(gg in unique(ts)){
    all_result[[gg]]=dm(gg,a,label)
    print(gg)
  }
  genus[a]=list(all_result)
  print(a)
}
names(genus)=c("method1","method2","method3")
rm(list=c("all_result"))

#########################################
##### under same family result ##########
#########################################
for(i in 1:length(family)){
  result=list()
  for(j in 1:length(table(data$family))){
    f_result=family[[i]][names(table(data$family))[j]]
    g_result=genus[[i]][tax_f_g(names(table(data$family))[j])]
    result[[j]]=list(f_result,g_result)
    names(result)[[j]]=names(table(data$family))[j]
    names(result[[j]])[1]='family'; names(result[[j]])[2]='genus'
    print(j)
  }
  #saveRDS(result,file=paste0(file,"result/",crop,"_",i,".rds"))
}
# rm(list=c("f_result","g_result","result"))
