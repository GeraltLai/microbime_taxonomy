#####
rm(list=ls())
#####
library(dada2)
packageVersion("dada2")
library(stringr)
# metadata
metadata=read.table(file="C:/Users/lab205/Desktop/wildrice_metadata.txt",sep=",",header=T)

##### parameter set set
file="C:/Users/lab205/Desktop/sra"
input=list.files(file,pattern = ".fastq$")
output_OTU="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/wildrice/ildrice_otu.txt"
output_taxa="C:/Users/lab205/Desktop/lab205/Microbiota_Analysis_BRC/wildrice/wildrice_taxa.txt"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs=sort(list.files(file, pattern="_1.fastq$",full.names = T))
fnRs=sort(list.files(file, pattern="_2.fastq$",full.names = T))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names=str_sub(input,1,nchar(input[1])-8)
sample.names=unique(sample.names)

### Inspect read quality profiles
# forward reads
plotQualityProfile(fnFs[1:2])
#  reverse reads
plotQualityProfile(fnRs[1:2])

### Filter and trim
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(file, "filtered", paste0(sample.names, "_1.fastq"))
filtRs <- file.path(file, "filtered", paste0(sample.names, "_2.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

### Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

### Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

### Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

### Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
nchar(mergers[[1]]$sequence[1])

### Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

### Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

### Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/lab205/Desktop/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/Users/lab205/Desktop/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

### Evaluate rownames(taxa) and colnames(seqtab.nochim)
for(i in 1:length(colnames(seqtab.nochim))){
  if(colnames(seqtab.nochim)[1]==rownames(taxa)[1]){
    print(i)
  }else{
    break
  }
}

### rename
name=lapply(1:length(colnames(seqtab.nochim)),function(x)paste0("OTU",x))
name=do.call(rbind,test)
colnames(seqtab.nochim)=name
rownames(taxa)=name
colnames(taxa)=c("kingdom", "phylum", "class", "order", "family", "genus", "species")
### save output
write.table(seqtab.nochim,file=output_OTU)
write.table(taxa,file=output_taxa)
test=read.table(file=output_taxa,header=T)
