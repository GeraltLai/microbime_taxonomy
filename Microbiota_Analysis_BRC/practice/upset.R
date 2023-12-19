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
all_result=readRDS(file=paste0(file,"result/",crop,"_group_3",".rds"))

#### pairwise reult
pairwise_result=readRDS(file=paste0(file,"result/",crop,"_pairwise.rds"))

#### minusone result
minusone_result=readRDS(file=paste0(file,"result/",crop,"_minusone.rds"))


####
i = 6
test1=as.data.frame(all_result[[i]])
test2=as.data.frame(pairwise_result[[i]])
test3=as.data.frame(minusone_result[[i]])
test=rbind(test1,test2,test3)
test=test[-1,]
genus=all_result[[i]][-1,]
genus_name=row.names(genus)
ggplot(test, aes(x=reorder(row.names(test),`Pr(>F)`), y= `Pr(>F)`)) +
  geom_bar(stat = "identity") +
  axis_combmatrix(sep = "-", levels = c(genus_name)) +
  geom_hline(aes(yintercept=0.05),linetype='dashed', colour = "blue")+
  xlab("Genus combination")+
  ylab("p-value")





##### test
test1=c("A","B","C","A-B","A-C","B-C","A-B-C")
test2=c(0.2,0.1,0.4,0.3,0.5,0.8,1.0)
test=as.data.frame(cbind(test1,test2))
test
ggplot(test, aes(x=test1, y=test2)) +
  geom_bar(stat = "identity")+
  geom_hline(aes(yintercept=0.05))+
  axis_combmatrix(sep = "-", levels = c("A","B","C")) 


###########################
###### ggupset plot #######
###########################
#total = rbind(step1_output, step2_output, step3_output)
#total = total[-1,]
#total
#genus = step1_output[-1,]
#genus_name = row.names(genus)
#ggplot(total, aes(x = reorder(row.names(total), `Pr(>F)`), y= `Pr(>F)`)) +
#  geom_bar(stat = "identity") +
#  axis_combmatrix(sep = "-", levels = c(genus_name)) +
#  geom_hline(aes(yintercept=0.05),linetype='dashed', colour = "blue")+
#  xlab("Genus combination")+
#  ylab("p-value")






