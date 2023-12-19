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
library(gtools)  # permutations packages
library(ggupset) # upset plot

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

#### Group3
summary1=readRDS(file=paste0(file,"result/",crop,"_group_3",".rds"))


#######################
###### function #######
#######################
pairwise=function(data,result,limit=0){
  com=list()
  for(a in 1:length(result)){
    all_output=c()
    name=c()
    pick=as.data.frame(result[[a]])
    pick1=subset(pick,pick[["size"]]>=max(pick[["size"]])*limit)
    for(i in 2:(nrow(pick1)-1))
      for(j in (i+1):nrow(pick1)){
        size=sum(pick1[i,1]+pick1[j,1])
        id1=rownames(pick1)[i]
        id2=rownames(pick1)[j]
        name=c(name,paste0(id1,"-",id2))
        new=data[data[["genus"]]==id1 | data[["genus"]]==id2,]
        clean=select(OTU,row.names(new))
        dist=bray0(clean)
        ans=adonis2(dist~meta[, label])[1,3:5]
        output=round(c(size,ans[1,1],ans[1,2],ans[1,3]),4)
        all_output=rbind(all_output,output)
      }
    colnames(all_output)=c("size","R2","F","Pr(>F)")
    row.names(all_output)=name
    com[a]=list(all_output)
    names(com)[a]=names(summary1[a])
  }
  return(com)
}
####
minus_one=function(data,result,limit=0){
  minus=list()
  for(a in 1:length(result)){
    all_output=c()
    name=c()
    pick=as.data.frame(result[[a]])
    id_f=rownames(pick)[1]
    new=data[data[["family"]]==id_f,]
    pick1=subset(pick,pick[["size"]]>=max(pick[["size"]])*limit)
    pick1=pick1[-1,]
    all_gene=row.names(pick1)
    for(i in row.names(pick1)){
      delet=which(row.names(pick1)==i)
      size=sum(pick1[-delet,1])
      part_gene=all_gene[-(delet)]
      name=c(name,paste0(part_gene,collapse="-"))
      new_g=new[new[["genus"]]!=i,]
      clean=select(OTU,row.names(new_g))
      dist=bray0(clean)
      ans=adonis2(dist~meta[, label])[1,3:5]
      output=round(c(size,ans[1,1],ans[1,2],ans[1,3]),4)
      all_output=rbind(all_output,output)
    }
    colnames(all_output)=c("size","R2","F","Pr(>F)")
    row.names(all_output)=name
    minus[a]=list(all_output)
    names(minus)[a]=names(result[a])
  }
  return(minus)
}


#### test
test=pairwise(data=data,result=summary1,limit=0)
# saveRDS(test,file=paste0(file,"result/",crop,"_pairwise.rds"))
test2=minus_one(data=data,result=summary1,limit=0)
# saveRDS(test2,file=paste0(file,"result/",crop,"_minusone.rds"))
