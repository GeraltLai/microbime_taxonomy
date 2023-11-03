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
library(cowplot) # combine plot

###### read data
rm(list = ls())
data = readRDS(file = "C:/Users/lab205/Desktop/basic.rds")
original = result_list
test = list()
for(i in 1:length(data)){
  test[[i]] = rbind(original[[i]], data[[i]])
}
names(test) = names(original)

data = test
saveRDS(data, file = "C:/Users/lab205/Desktop/test/basic.rds")

###### read data
rm(list = ls())
data = readRDS(file = "C:/Users/lab205/Desktop/method4_case_less_than_control.rds")
data = readRDS(file = "C:/Users/lab205/Desktop/method4_case_more_than_control.rds")

### genus1
genus1 = data[[2]]
### significant result
data1 = genus1[which(genus1$individual == 1),]
# total
test1 = data1[,c("total","delta1")]
p1 = ggplot(data = test1, aes(x = as.factor(delta1) , fill = as.factor(total)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "significant result", values = c("#E69F00", "#56B4E9", "#009E73", 
                                                            "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage significant result") +
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("significant genus(TP)") +
  theme(plot.title = element_text(hjust = 0.5))
p1

# with_sig
test2 = data1[,c("with_sig","delta1")]
p2 = ggplot(data = test2, aes(x = as.factor(delta1) , fill = as.factor(with_sig)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "with_sig_g", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                                    "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage with significant genus") +
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("significant genus(TP)") +
  theme(plot.title = element_text(hjust = 0.5))
p2

# with_non_sig
test3 = data1[,c("with_nonsig","delta1")]
p3 = ggplot(data = test3, aes(x = as.factor(delta1) , fill = as.factor(with_nonsig)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "with_nonsig_g", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                                       "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage with nonsignificant genus") +
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("significant genus(TP)") +
  theme(plot.title = element_text(hjust = 0.5))
p3

### nonsignificant result 
data2 = genus1[which(genus1$individual == 0),]
# total
test4 = data2[,c("total","delta1")]
p4 = ggplot(data = test4, aes(x = as.factor(delta1) , fill = as.factor(total)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "significant result", values = c("#999999","#E69F00", "#56B4E9", "#009E73", 
                                                            "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage significant result") +
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("nonsignificant genus(FN)") +
  theme(plot.title = element_text(hjust = 0.5))
p4

# with_sig
test5 = data2[,c("with_sig","delta1")]
p5 = ggplot(data = test5, aes(x = as.factor(delta1) , fill = as.factor(with_sig)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "with_sig_g", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                                    "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage with significant genus") +
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("nonsignificant genus(FN)") +
  theme(plot.title = element_text(hjust = 0.5))
p5

# with_non_sig
test6 = data2[,c("with_nonsig","delta1")]
p6 = ggplot(data = test6, aes(x = as.factor(delta1) , fill = as.factor(with_nonsig)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "with_nonsig_g", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                                       "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage with nonsignificant genus") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("nonsignificant genus(FN)") +
  theme(plot.title = element_text(hjust = 0.5))
p6

### genus4
genus4 = data[[7]]
### significant result
data3 = genus4[which(genus4$individual == 1),]
# total
test7 = data3[,c("total","delta1")]
p7 = ggplot(data = test7, aes(x = as.factor(delta1) , fill = as.factor(total)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "significant result", values = c("#E69F00", "#56B4E9", "#009E73", 
                                                            "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage significant result") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("significant genus(FP)") +
  theme(plot.title = element_text(hjust = 0.5))
p7

# with_sig
test8 = data3[,c("with_sig","delta1")]
p8 = ggplot(data = test8, aes(x = as.factor(delta1) , fill = as.factor(with_sig)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "with_sig_g", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                                    "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x = "delta1", y = "Percentage with significant genus") +
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("significant genus(FP)") +
  theme(plot.title = element_text(hjust = 0.5))
p8

# with_nonsig
test9 = data3[,c("with_nonsig","delta1")]
p9 = ggplot(data = test9, aes(x = as.factor(delta1) , fill = as.factor(with_nonsig)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "with_nonsig_g", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                                       "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x = "delta1", y = "Percentage with nonsignificant genus") +
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("significant genus(FP)") +
  theme(plot.title = element_text(hjust = 0.5))
p9

### nonsignificant result 
data4 = genus4[which(genus4$individual == 0),]
# total
test10 = data4[,c("total","delta1")]
p10 = ggplot(data = test10, aes(x = as.factor(delta1) , fill = as.factor(total)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "significant result", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                                            "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage significant result") +
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("nonsignificant genus(TN)") +
  theme(plot.title = element_text(hjust = 0.5))
p10

# with_sig
test11 = data4[,c("with_sig","delta1")]
p11 = ggplot(data = test11, aes(x = as.factor(delta1) , fill = as.factor(with_sig)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "with_sig_g", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                                    "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage with significant genus") +
  scale_y_continuous(labels = scales::percent) + 
  ggtitle("nonsignificant genus(TN)") +
  theme(plot.title = element_text(hjust = 0.5))
p11

# with_non_sig
test12 = data4[,c("with_nonsig","delta1")]
p12 = ggplot(data = test12, aes(x = as.factor(delta1) , fill = as.factor(with_nonsig)))+
  geom_bar(stat="count", position ="fill") +
  scale_fill_manual(name = "with_nonsig_g", values = c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                                                       "#F0E442", "#0072B2", "darkorchid1","#D55E00", "#CC79A7","red"))+
  labs(x ="delta1", y = "Percentage with nonsignificant genus") +
  scale_y_continuous(labels = scales::percent) +
  ggtitle("nonsignificant genus(TN)") +
  theme(plot.title = element_text(hjust = 0.5))
p12


### case : control = 1 :2
TP_FP_total =plot_grid(p1, p7, labels = "AUTO")
TP_FP_with_sig =plot_grid(p2, p8, labels = "AUTO")
TP_FP_with_nonsig =plot_grid(p3, p9, labels = "AUTO")
TN_FN_total =plot_grid(p4, p10, labels = "AUTO")
TN_FN_wtih_sig = plot_grid(p5, p11, labels = "AUTO")
TN_FN_wtih_nonsig = plot_grid(p6, p12, labels = "AUTO")

### case : control = 1 :2
title <- ggdraw() + draw_label("case studies : control studies = 1:2" , size = 20, fontface='bold')
plot1 = plot_grid(title, TP_FP_total, ncol=1, rel_heights=c(0.1, 1))
plot2 = plot_grid(title, TN_FN_total, ncol=1, rel_heights=c(0.1, 1))

### case : control = 2 :1
TP_FP_total = plot_grid(p1, p7, labels = c("C", "D"))
TP_FP_with_sig = plot_grid(p2, p8, labels = c("C", "D"))
TP_FP_with_nonsig = plot_grid(p3, p9, labels = c("C", "D"))
TN_FN_total = plot_grid(p4, p10, labels = c("C", "D"))
TN_FN_wtih_sig = plot_grid(p5, p11, labels = c("C", "D"))
TN_FN_wtih_nonsig = plot_grid(p6, p12, labels = c("C", "D"))

### case : control = 2 :1
title <- ggdraw() + draw_label("case studies : control studies = 2:1" , size = 20, fontface='bold')
plot3 = plot_grid(title, TP_FP_total, ncol=1, rel_heights=c(0.1, 1))
plot4 = plot_grid(title, TN_FN_total, ncol=1, rel_heights=c(0.1, 1))

### combine different scenario
pp1 = plot_grid(plot1, plot3, labels = NULL, ncol = 1)
pp2 = plot_grid(plot2, plot4, labels = NULL, ncol = 1)
pp1
pp2

### case = control
pp1 = plot_grid(p1, p7, p4, p10, labels = "AUTO", ncol = 2)
pp1
pp2 = plot_grid(p2, p8, p5, p11, labels = "AUTO", ncol = 2)
pp2
pp3 = plot_grid(p3, p9, p6, p12, labels = "AUTO", ncol = 2)
pp3


# save plot
ggsave(
  filename = "method4.png",
  plot = pp1,
  path = "C:/Users/lab205/Desktop/test",
  scale = 1,
  width = NA,
  height = NA,
  dpi = 300,
  bg = "white"
)
"#E69F00"




###
title <- ggdraw() + draw_label("case studies : control studies = 1:2" , size = 20, fontface='bold')
plot1= plot_grid(title, TP_FP_total, ncol=1, rel_heights=c(0.1, 1))
plot2 = plot_grid(title, TN_FN_total, ncol=1, rel_heights=c(0.1, 1))

###
title <- ggdraw() + draw_label("case studies : control studies = 2:1" , size = 20, fontface='bold')
plot3 = plot_grid(title, TP_FP_total, ncol=1, rel_heights=c(0.1, 1))
plot4 = plot_grid(title, TN_FN_total, ncol=1, rel_heights=c(0.1, 1))




