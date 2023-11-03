#####
rm(list=ls())
library(stringr) 

#####
file="C:/Users/lab205/Desktop/sra"
input=list.files(file,pattern = ".fastq$")
input

####
sample=str_sub(input,1,nchar(input[1])-8) 
direction=rep(c("forward","reverse"),length(input)/2)
absolute=paste0("$PWD/",input)
data=data.frame(sample,absolute,direction)
names(data)=c('sample-id','absolute-filepath','direction')
write.csv(x=data,file="C:/Users/lab205/Desktop/manifest.csv",row.names=F,quote=F)
