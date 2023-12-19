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

#####################
##### load data #####
#####################
OTU = read.table("C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/tutorial/example.final.an.unique_list.0.03.norm.shared.txt", header=TRUE, sep="\t")

# Taxonomy of each OTU
tax = read.table("C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/tutorial/example.final.an.unique_list.0.03.cons.taxonomy.txt", header=TRUE, sep="\t")

# Metadata. Since we made this in Excel, not mothur, we can use the "row.names" modifier to automatically name the rows by the values in the first column (sample names)
meta = read.table("C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/tutorial/example.metadata.txt", header=TRUE, row.names=1, sep="\t")

# SCFA data
SCFA = read.table("C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/tutorial/example.SCFA.txt", header=TRUE, row.names=1, sep="\t")

#############################
##### Clean up the data #####
#############################
# OTU table
row.names(OTU) = OTU$Group
OTU.clean = OTU[,-which(names(OTU) %in% c("label", "numOtus", "Group"))]

# Taxonomy table
row.names(tax) = tax$OTU
tax.clean = tax[row.names(tax) %in% colnames(OTU.clean),]
tax.clean = separate(tax.clean, Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"), sep=";")
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("Size", "Strain", "OTU"))]

# Order the data
OTU.clean = OTU.clean[order(row.names(OTU.clean)),]
meta = meta[order(row.names(meta)),]
SCFA = SCFA[order(row.names(SCFA)),]

# set seed
set.seed(8765)

###########################
##### Alpha-diversity #####
###########################
### Explore alpha metrics
par(mfrow = c(2, 2))
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(meta$simpson, main="Simpson diversity", xlab="", breaks=10)
hist(meta$chao, main="Chao richness", xlab="", breaks=15)
hist(meta$ace, main="ACE richness", xlab="", breaks=15)

# Create 2x2 plot environment 
par(mfrow = c(2, 2))
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(1/meta$simpson, main="Inverse Simpson diversity", xlab="", breaks=10)
hist(meta$chao, main="Chao richness", xlab="", breaks=15)
hist(meta$ace, main="ACE richness", xlab="", breaks=15)
# Shapiro-Wilk normality test
shapiro.test(meta$shannon)
shapiro.test(1/meta$simpson)
shapiro.test(meta$chao)
shapiro.test(meta$ace)

### Categorical variables
aov.shannon.age = aov(shannon ~ AgeGroup, data=meta)
summary(aov.shannon.age)
# Tukey HSD
TukeyHSD(aov.shannon.age)
# Re-order the groups because the default is 1yr-2w-8w
meta$AgeGroup.ord = factor(meta$AgeGroup, c("2w","8w","1yr"))
# Return the plot area to 1x1
par(mfrow = c(1, 1))
boxplot(shannon ~ AgeGroup.ord, data=meta, ylab="Shannon's diversity")

# Non-normally distributed metrics
kruskal.test(chao ~ AgeGroup, data=meta)
pairwise.wilcox.test(meta$chao, meta$AgeGroup, p.adjust.method="fdr")
# Plot
par(mfrow = c(1, 1))
boxplot(chao ~ AgeGroup.ord, data=meta, ylab="Chao richness")

### Continuous variables
# Normally distributed metrics
glm.shannon.ADG = glm(shannon ~ ADGKG, data=meta)
summary(glm.shannon.ADG)
plot(shannon ~ ADGKG, data=meta)
abline(glm.shannon.ADG)

# Non-normally distributed metrics
gaussian.chao.ADG = glm(chao ~ ADGKG, data=meta, family="gaussian")
par(mfrow = c(1,2))
plot(gaussian.chao.ADG, which=c(1,2))
qp.chao.ADG = glm(chao ~ ADGKG, data=meta, family="quasipoisson")
par(mfrow = c(1,2))
plot(qp.chao.ADG, which=c(1,2))
summary(qp.chao.ADG)

# Plot
par(mfrow = c(1, 1))
plot(log(chao) ~ ADGKG, data=meta, ylab="ln(Chao's richness)")
abline(qp.chao.ADG)

### Mixed models
aov.shannon.all = aov(shannon ~ AgeGroup*ADGKG, data=meta)
summary(aov.shannon.all)
aov.shannon.all2 = aov(shannon ~ AgeGroup+ADGKG, data=meta)
summary(aov.shannon.all2)
TukeyHSD(aov.shannon.all)
glm.shannon.all = glm(shannon ~ AgeGroup*ADGKG, data=meta)
summary(glm.shannon.all)
glm.shannon.all2 = glm(shannon ~ AgeGroup+ADGKG, data=meta)
summary(glm.shannon.all2)
qp.chao.all = glm(chao ~ AgeGroup*ADGKG, data=meta, family="quasipoisson")
summary(qp.chao.all)
qp.chao.all2 = glm(chao ~ AgeGroup+ADGKG, data=meta, family="quasipoisson")
summary(qp.chao.all2)

### Repeated measure
rm.shannon.all = lmer(shannon ~ AgeGroup+ADGKG + (1|Animal), data=meta)
summary(rm.shannon.all)

##########################
##### Beta-diversity #####
##########################
### OTU-based metrics
BC.nmds = metaMDS(OTU.clean, distance="bray", k=2, trymax=1000)
# plot
par(mfrow = c(1, 1))
plot(BC.nmds, type="n", main="Bray-Curtis",xlim=c(-6,4),ylim=c(-3,2))
points(BC.nmds, display="sites", pch=20, col=c("green","red","blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])
legend(-5.5, 2.5, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)

# Jaccard metric
J.nmds = metaMDS(OTU.clean, distance="jaccard", k=2, trymax=1000)
plot(J.nmds, type="n", main="Jaccard",xlim=c(-3,2),ylim=c(-1.5,1.5))
points(J.nmds, display="sites", pch=20, col=c("green","red","blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])
legend(-3, 1.5, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)

### Ellipses
plot(BC.nmds, type="n", main="Bray-Curtis",xlim=c(-6,4),ylim=c(-3,2))
legend(-5.5, 3, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)
#Add an ellipse for 2w
ordiellipse(BC.nmds, groups=factor(meta$AgeGroup, c("2w","8w","1yr")), display="sites", kind="se", conf=0.99, label=FALSE, col="green", draw="polygon", alpha=200, show.groups = c("2w"), border=FALSE)
#Add an ellipse for 8w
ordiellipse(BC.nmds, groups=factor(meta$AgeGroup, c("2w","8w","1yr")), display="sites", kind="se", conf=0.99, label=FALSE, col="red", draw="polygon", alpha=200, show.groups = c("8w"), border=FALSE)
#Add an ellipse for 1yr
ordiellipse(BC.nmds, groups=factor(meta$AgeGroup, c("2w","8w","1yr")), display="sites", kind="se", conf=0.99, label=FALSE, col="blue", draw="polygon", alpha=200, show.groups = c("1yr"), border=FALSE)

### Extract x-y-z values for this nmds
BC.nmds.3D = metaMDS(OTU.clean, distance="bray", k=3, trymax=1000)
BCxyz = scores(BC.nmds.3D, display="sites")
BCxyz
# 3D plot
plot_ly(x=BCxyz[,1], y=BCxyz[,2], z=BCxyz[,3], type="scatter3d", mode="markers", color=factor(meta$AgeGroup, c("2w","8w","1yr")), colors=c("green", "red", "blue"))
# plot
par(mfrow=c(1,2))
# Axis 1 and 2 (x and y)
plot(BCxyz[,1], BCxyz[,2], main="Bray-Curtis 1:2", pch=20, col=c("green","red","blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])
legend(-5, 3, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)
#Axis 1 and 3 (x and z)
plot(BCxyz[,1], BCxyz[,3], main="Bray-Curtis 1:3", pch=20, col=c("green","red","blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])

### Phylogentic-based metrics
# Create physeq object
OTU.UF = otu_table(as.matrix(OTU.clean), taxa_are_rows=FALSE)
tax.UF = tax_table(as.matrix(tax.clean))
meta.UF = sample_data(meta)
# merge
physeq = phyloseq(OTU.UF, tax.UF, meta.UF)
# load 
load("C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/NJ.tree.Rdata")
# UniFrac calculations
physeq.tree = merge_phyloseq(physeq, NJ.tree)
physeq.tree
# Dot plots
wUF.ordu = ordinate(physeq.tree, method="NMDS", distance="unifrac", weighted=TRUE)
par(mfrow=c(1,1))
plot(wUF.ordu, type="n", main="Weighted UniFrac")
points(wUF.ordu, pch=20, display="sites", col=c("green","red","blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])
legend(0.3,0.15, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)
# ggplot2 package
plot_ordination(physeq.tree, wUF.ordu, type="sites", color="AgeGroup") + 
  scale_colour_manual(values=c("2w"="green", "8w"="red", "1yr"="blue")) + 
  theme_bw() + 
  ggtitle("Weighted UniFrac")

# Unweighted UniFrac 
uwUF.ordu = ordinate(physeq.tree, method="NMDS", distance="unifrac", weighted=FALSE)
plot_ordination(physeq.tree, uwUF.ordu, type="sites", color="AgeGroup") + 
  scale_colour_manual(values=c("2w"="green", "8w"="red", "1yr"="blue")) + 
  theme_bw() + 
  ggtitle("Unweighted UniFrac")

# Ellipses
plot(wUF.ordu, type="n", main="Weighted UniFrac")
legend(0.3, 0.15, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)

# Add an ellipse for 2w
ordiellipse(wUF.ordu, groups=factor(meta$AgeGroup, c("2w","8w","1yr")), display="sites", kind="se", conf=0.99, label=FALSE, col="green", draw="polygon", alpha=200, show.groups = c("2w"), border=FALSE)

# Add an ellipse for 8w
ordiellipse(wUF.ordu, groups=factor(meta$AgeGroup, c("2w","8w","1yr")), display="sites", kind="se", conf=0.99, label=FALSE, col="red", draw="polygon", alpha=200, show.groups = c("8w"), border=FALSE)

# Add an ellipse for 1yr
ordiellipse(wUF.ordu, groups=factor(meta$AgeGroup, c("2w","8w","1yr")), display="sites", kind="se", conf=0.99, label=FALSE, col="blue", draw="polygon", alpha=200, show.groups = c("1yr"), border=FALSE)

# plot ellipses with ggplot2 by adding the stat_ellipse
plot_ordination(physeq.tree, wUF.ordu, type="sites", color="AgeGroup") + 
  scale_colour_manual(values=c("2w"="green", "8w"="red", "1yr"="blue")) + 
  theme_bw() + 
  stat_ellipse() + 
  ggtitle("Weighted UniFrac")

# 3D plots
wUF.ordu
wUF.ordu = UniFrac(physeq.tree, weighted=TRUE, normalized=TRUE)
wUF.nmds.3D = metaMDS(wUF.ordu, method="NMDS", k=3)
wUFxyz = scores(wUF.nmds.3D, display="sites")
wUFxyz
plot_ly(x=wUFxyz[,1], y=wUFxyz[,2], z=wUFxyz[,3], type="scatter3d", mode="markers", color=factor(meta$AgeGroup, c("2w","8w","1yr")), colors=c("green", "red", "blue"))

##### Vectors for continuous variables #####
fit.BC = envfit(BC.nmds, meta) 
fit.BC
fit.BC = envfit(BC.nmds, meta[,c("AgeGroup", "ADGKG")])
fit.BC
fit.wUF = envfit(wUF.ordu, meta[,c("AgeGroup", "ADGKG")])
fit.wUF

# plot
plot(BC.nmds, type="n", main="Bray-Curtis",xlim=c(-6,4),ylim=c(-3,2))
points(BC.nmds, pch=20, display="sites", col=c("green", "red", "blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])
legend(-6, 2, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)
# Add fitted variables
plot(fit.BC, col="black")

# only plot variables with a fit P-value < 0.05
plot(BC.nmds, type="n", main="Bray-Curtis",xlim=c(-6,4),ylim=c(-3,2))
points(BC.nmds, pch=20, display="sites", col=c("green", "red", "blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])
legend(-6, 2, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)
# Add fitted variables
plot(fit.BC, col="black", p.max=0.05)

# Weighted UniFrac
plot(wUF.ordu, type="n", main="Weighted UniFrac")
points(wUF.ordu, pch=20, display="sites", col=c("green", "red", "blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])
legend(.3,.15, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)
# Add fitted variables
plot(fit.wUF, col="black")

# Fitting all OTUs would take awhile so we will only fit the first 10 in our table.
fit.BC.OTU = envfit(BC.nmds, OTU.clean[,1:10])
fit.BC.OTU
# We will only plot significant arrows in this case
plot(BC.nmds, type="n", main="Bray-Curtis",xlim=c(-6,4),ylim=c(-3,2))
points(BC.nmds, pch=20, display="sites", col=c("green", "red", "blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])
legend(-6, -1.1, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)
# Add fitted variables
plot(fit.BC.OTU, col="black", p.max=0.05)

# Extract all OTUs within the genus Ruminococcus
OTU.Rumino = OTU.clean[,tax.clean$Genus == "g__Ruminococcus"]
# Sum the abundances of the Ruminococcaceae OTUs into one variable (column)
OTU.Rumino$Rumino.sum = rowSums(OTU.Rumino)
# Fit the new Ruminococcaceae group
fit.BC.Rumino = envfit(BC.nmds, OTU.Rumino$Rumino.sum)
fit.BC.Rumino

# Plot
plot(BC.nmds, type="n", main="Bray-Curtis",xlim=c(-6,4),ylim=c(-3,2))
points(BC.nmds, pch=20, display="sites", col=c("green", "red", "blue")[factor(meta$AgeGroup, c("2w","8w","1yr"))])
legend(-6, -1.1, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)
# Add fitted variables
plot(fit.BC.Rumino, col="black", labels=c("Ruminococcus"))

##### PERMANOVA #####
# Calculate distance and save as a matrix
BC.dist=vegdist(OTU.clean, distance="bray")
# Run PERMANOVA on distances.
adonis(BC.dist ~ AgeGroup*ADGKG, data = meta, permutations = 1000)
# Similarly for Jaccard
J.dist=vegdist(OTU.clean, distance="jaccard")
adonis(J.dist ~ AgeGroup*ADGKG, data = meta, permutations = 1000)
# We see that the interaction is not significant so we remove it.
adonis(BC.dist ~ AgeGroup+ADGKG, data = meta, permutations = 1000)
adonis(J.dist ~ AgeGroup+ADGKG, data = meta, permutations = 1000)

# We use the phyloseq package to calculate distances and then vegan to run PERMANOVA.
wUF.dist = UniFrac(physeq.tree, weighted=TRUE, normalized=TRUE)
adonis(wUF.dist ~ AgeGroup*ADGKG, data=meta, permutations = 1000)
uwUF.dist = UniFrac(physeq.tree, weighted=FALSE, normalized=TRUE)
adonis(uwUF.dist ~ AgeGroup*ADGKG, data=meta, permutations = 1000)
# Remove non-significant interaction term
adonis(wUF.dist ~ AgeGroup+ADGKG, data=meta, permutations = 1000)
adonis(uwUF.dist ~ AgeGroup+ADGKG, data=meta, permutations = 1000)

##### ANOSIM #####
# Bray-Curtis
anosim(BC.dist, meta$AgeGroup, permutations = 1000)

##### 2D variables #####
OTU.SCFA = OTU.clean[row.names(OTU.clean) %in% paste(row.names(SCFA), ".F", sep=""),]
dist1 = vegdist(OTU.SCFA)
dist2 = vegdist(SCFA)
mantel(dist1, dist2, permutations=100)
# Run a Mantel test comparing the 2 matrices.
mantel(dist1, dist2, permutations=100)

##### Beta dispersion #####
disp.age = betadisper(BC.dist, meta$AgeGroup)
permutest(disp.age, pairwise=TRUE, permutations=1000)
# plot
plot(BC.nmds, type="n", main="Bray-Curtis",xlim=c(-6,4),ylim=c(-3,2))
legend(.6,-2, legend=c("2w","8w","1yr"), col=c("green","red","blue"), pch=20)
ordiellipse(BC.nmds, groups=factor(meta$AgeGroup, c("2w","8w","1yr")), display="sites", kind="se", conf=0.99, label=FALSE, col="green", draw="polygon", alpha=200, show.groups = c("2w"), border=FALSE)
ordiellipse(BC.nmds, groups=factor(meta$AgeGroup, c("2w","8w","1yr")), display="sites", kind="se", conf=0.99, label=FALSE, col="red", draw="polygon", alpha=200, show.groups = c("8w"), border=FALSE)
ordiellipse(BC.nmds, groups=factor(meta$AgeGroup, c("2w","8w","1yr")), display="sites", kind="se", conf=0.99, label=FALSE, col="blue", draw="polygon", alpha=200, show.groups = c("1yr"), border=FALSE)


