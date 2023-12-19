#####
rm(list=ls())

#####  packages & function
library(vegan)
library(testthat)
source(file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/permanova.R")

####################
##### dune data ####
####################
data(dune)
data(dune.env)
data1=dune
meta1=dune.env

# sum of square test
mg=manova(data=data1,meta=meta1,factor="Management")
dune.dist <- vegdist(data1, method="bray")
pg <- adonis2(dune.dist ~ Management, data = meta1, permutations = 999, method="bray")
test_that("manova",{
  expect_equal(as.numeric(c(mg[1],mg[2],mg[3],mg[4])),
               c(pg[1,2],pg[2,2],pg[3,2],pg[1,4]))
})
# p-value test
myf=permutation(f_value=mg[4],data=data1,meta=meta1,factor="Management")
test_that("permanova of p_value",{
  expect_equal(as.numeric(myf),pg[1,5])
})


################
### cow data ###
################
data3=read.table(file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/cow/cow_otu.txt",header=T)
meta3=read.table(file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/cow/cow_metadata.txt",header=T)

# sum of square test
mg=manova(data=data3,meta=meta3,factor="AgeGroup")
dune.dist <- vegdist(data3, method="bray")
pg <- adonis2(dune.dist ~ AgeGroup, data = meta3, permutations = 999, method="bray")
test_that("manova",{
  expect_equal(as.numeric(c(mg[1],mg[2],mg[3],mg[4])),
               c(pg[1,2],pg[2,2],pg[3,2],pg[1,4]))
})
# p-value test
myf=permutation(f_value=mg[4],data=data3,meta=meta3,factor="AgeGroup")
test_that("permanova of p_value",{
  expect_equal(as.numeric(myf),pg[1,5])
})


#######################
### wildrice data #####
#######################
data2=read.table(file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/wildrice/wildrice_otu.txt",header=T)
meta2=read.table(file="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/wildrice/wildrice_metadata.txt",header=T)

# sum of square test
mg=manova(data=data2,meta=meta2,factor="experimental_factor")
dune.dist <- vegdist(data2, method="bray")
pg <- adonis2(dune.dist ~ experimental_factor, data = meta2, permutations = 999, method="bray")
test_that("manova",{
  expect_equal(as.numeric(c(mg[1],mg[2],mg[3],mg[4])),
               c(pg[1,2],pg[2,2],pg[3,2],pg[1,4]))
})
# p-value test
myf=permutation(f_value=mg[4],data=data2,meta=meta2,factor="experimental_factor")
test_that("permanova of p_value",{
  expect_equal(as.numeric(myf),pg[1,5])
})



