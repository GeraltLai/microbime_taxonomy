#####
rm(list=ls())
library(UpSetR)
library(ComplexHeatmap)
library(pheatmap)
library(stringr)
#####
result1=readRDS(file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/wildrice/result/wildrice_group_3.rds")
result2=readRDS(file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/wildrice/result/wildrice_combine_2.rds")
result3=readRDS(file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/wildrice/result/wildrice_minus_1.rds")

#####
i=names(result1)[3]
result1[i]
result2[i]
result3[i]

#####
asdfg <- names(result1)[4]
result1[i]
result2[i]
result3[i]
# Incertae Sedis p-value<0.05
# Lachnotalea p-value>0.05
# Incertae Sedis>Lachnotalea

#####
i=names(result1)[5]
result1[i]
result2[i]
result3[i]

#####
i=names(result1)[6]
result1[i]
result2[i]
result3[i]

### heatmap
test1=result1[[4]]
test1=test1[-1,]
test2=result2[[4]]
test3=result3[[4]]

###
num_col <- nrow(test1) + nrow(test2) + nrow(test3)
m <- matrix(data=0, nrow=nrow(test1), ncol=num_col)

name=c(row.names(test1),row.names(test2),row.names(test3))
row.names(m)=row.names(test1)
colnames(m)=name
for(i in 1:nrow(m)){
  m[i,i]=test1[i,4]
}

for(i in 1:nrow(test2)){
  n1=row.names(test2)[i]
  n2=str_split(n1,pattern = "_&_")
  n3=n2[[1]][1]
  n4=n2[[1]][2]
  m[n3,i+nrow(test1)]=test2[i,4]
  m[n4,i+nrow(test1)]=test2[i,4]
}

for(i in 1:nrow(test3)){
  n1=row.names(test3)[i]
  n2=str_split(n1,pattern = "minus_")
  n3=n2[[1]][2]
  n4=which(row.names(m)==n3)
  m[-n4,i+nrow(test1)+nrow(test2)]=test3[i,4]
  
}
### default value
pheatmap(m)

### different color
pheatmap(m, color = colorRampPalette(c("navy", "white", "firebrick3"))(5))

### tag the value
pheatmap(m, display_numbers = TRUE)






