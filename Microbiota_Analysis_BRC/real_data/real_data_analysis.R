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

### method1
m1=readRDS(file=paste0(file,"result/",crop,"_1.rds"))
### method2
m2=readRDS(file=paste0(file,"result/",crop,"_2.rds"))
### method3
m3=readRDS(file=paste0(file,"result/",crop,"_3.rds"))
### data
data=read.table(file=paste0(file,crop,"_data.txt"))

####################
##### function #####
####################
tax_result=function(x, fs=c(T, F), gs=c(T, F), alpha=0.05){
  result1=list();result2=list();result3=list();result4=list()
  for(i in 1:length(x)){
    print(i)
    f_result=((x[[i]][[1]][1])[[1]][1,3])
    g_result=x[[i]][[2]]
    count=0
    if(f_result<=alpha){
      for(j in 1:length(g_result)){
        if(g_result[[j]][1,3]<=alpha){
          count=count+1
        }else{
          count=count
        }
      }
      if(count==length(g_result)){
        result1=c(result1,x[i])
      }else{
        result3=c(result3,x[i])
      }
    }
    if(f_result>alpha){
      for(j in 1:length(g_result)){
        if(g_result[[j]][1,3]>alpha){
          count<-count+1
        }else{
          count<-count
        }
      }
      if(count==length(g_result)){
        result2=c(result2,x[i])
      }else{
        result4=c(result4,x[i])
      }
    }
  }
  result=list(result1,result2,result3,result4)
  names(result)=c("fs=T,gs=T","fs=F,gs=F","fs=T,gs=F","fs=F,gs=T")
  if(length(fs)+length(ts)>2) return(result)
  if(fs==T && gs==T) return(result1)
  if(fs==F && gs==F) return(result2)
  if(fs==T && gs==F) return(result3)
  if(fs==F && gs==T) return(result4)
}

###########################
##### analysis result #####
###########################
### Na value need notice
output=tax_result(m3)
output1=tax_result(m3,fs=T,gs=T,alpha=0.05)
output2=tax_result(m3,fs=F,gs=F,alpha=0.05)
output3=tax_result(m3,fs=T,gs=F,alpha=0.05)
output4=tax_result(m3,fs=F,gs=T,alpha=0.05)
sum=length(output1)+length(output2)+length(output3)+length(output4)
sum


################################################################
##### family significant and genus no significant analysis #####
################################################################
summary1=list()
for(i in 1:length(output3)){
  print(i)
  result=c()
  f_result=(output3[[i]][[1]])
  f_data=subset(data,data$family==names(f_result))
  f_size=sum(f_data$Size)
  result=rbind(result,c(f_size,f_result[[1]][1,1],f_result[[1]][1,2],f_result[[1]][1,3]))
  g_result=(output3[[i]][[2]])
  for(j in 1:length(g_result)){
    g_data=subset(f_data,f_data$genus==names(g_result)[j])
    g_size=sum(g_data$Size)
    result=rbind(result,c(g_size,g_result[[j]][1,1],g_result[[j]][1,2],g_result[[j]][1,3]))
  }
  result=round(result,4)
  row.names(result)=c(names(f_result),names(g_result))
  colnames(result)=c("size","R2","F","Pr(>F)")
  summary1[[i]]=result
  names(summary1)[[i]]=names(f_result)
}
summary1


################################################################
##### family no significant and genus significant analysis #####
################################################################
summary2=list()
for(i in 1:length(output4)){
  print(i)
  result=c()
  f_result=(output4[[i]][[1]])
  f_data=subset(data,data$family==names(f_result))
  f_size=sum(f_data$Size)
  result=rbind(result,c(f_size,f_result[[1]][1,1],f_result[[1]][1,2],f_result[[1]][1,3]))
  g_result=(output4[[i]][[2]])
  for(j in 1:length(g_result)){
    g_data=subset(f_data,f_data$genus==names(g_result)[j])
    g_size=sum(g_data$Size)
    result=rbind(result,c(g_size,g_result[[j]][1,1],g_result[[j]][1,2],g_result[[j]][1,3]))
  }
  result=round(result,4)
  row.names(result)=c(names(f_result),names(g_result))
  colnames(result)=c("size","R2","F","Pr(>F)")
  summary2[[i]]=result
  names(summary2)[[i]]=names(f_result)
}
summary2


###############################################
##### pick up summary1 interesting result #####
###############################################
id=c()
for(i in 1:length(summary1)){
  test=summary1[[i]]
  test1=test[order(-test[,1]),]
  if(test1[1,ncol(test1)]<=0.05 & test1[2,ncol(test1)]<=0.05){
    id=c(id,i)
  }
}
int=summary1[-id]
for(i in 1:length(int)){
  int[[i]]=int[[i]][order(-int[[i]][,1]),]
}
int

# saveRDS(int,file=paste0(file,"result/",crop,"_group_3",".rds"))


