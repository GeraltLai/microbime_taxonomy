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
meta = read.table(paste0(file,"wildrice_metadata.txt"), header=TRUE, sep=",")

### according to different part of plants
ex1=meta[meta[["experimental_factor"]]=="s",]
ex2=meta[meta[["experimental_factor"]]=="l",]
ex3=meta[meta[["experimental_factor"]]=="r",]
ex4=meta[meta[["experimental_factor"]]=="rh",]
ex5=meta[meta[["experimental_factor"]]=="t",]
id1=which(meta[["experimental_factor"]]=="s")
id2=which(meta[["experimental_factor"]]=="l")
OTU=OTU[c(id1,id2),]
meta=rbind(ex1,ex2)
size=lapply(1:ncol(OTU),function(x) sum(OTU[,x]))
size=do.call(rbind,size)
data=cbind(size,tax)
data=data[,-ncol(data)]
data=data[complete.cases(data),]


#######################
###### function #######
#######################
### Bray–Curtis distance matrix 
dm=function(x, method){
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
    return(adonis2(dist~Altitude*AvgSpotLen,data=meta)[1:3,3:5])
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
    if((length(table(new_meta$Altitude))<=1 | length(table(new_meta$AvgSpotLen))<=1)){
      return(matrix(rep(0,9),nrow=3,ncol=3))
    }else{
      dist=vegdist(clean, distance="bray")
      return(adonis2(dist~Altitude*AvgSpotLen,data=new_meta)[1:3,3:5])
    }
  }
  if(method==3){
    dist=bray0(clean)
    return(adonis2(dist~Altitude*AvgSpotLen,data=meta)[1:3,3:5])
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

family=list()
ts=data$family
new=data[data[["family"]]==unique(ts)[3],]
clean=select(OTU,row.names(new))
dist=bray0(clean)


#############################
##### family result ##########
#############################
family=list()
ts=data$family
for(a in 3:3){
  all_result=list()
  for(ff in unique(ts)){
    all_result[[ff]]=dm(ff,a)
    print(ff)
  }
  family[a]=list(all_result)
  print(a)
}
names(family)="method3"
rm(list=c("all_result"))
# saveRDS(family,file=paste0(file,"analysis/family.rds"))

#############################
##### genus result ##########
#############################
genus=list()
ts=data$genus
for(a in 1:3){
  all_result=list()
  for(gg in unique(ts)){
    all_result[[gg]]=dm(gg,a)
    print(gg)
  }
  genus[a]=list(all_result)
  print(a)
}
names(genus)=c("method1","method2","method3")
rm(list=c("all_result"))
# saveRDS(genus,file=paste0(file,"analysis/genus.rds"))

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
  saveRDS(result,file=paste0(file,"result/wildrice_",i,".rds"))
}
#rm(list=c("f_result","g_result","result"))




