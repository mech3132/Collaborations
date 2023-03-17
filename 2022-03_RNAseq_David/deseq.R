#!bin/bash Rscript
library("DESeq2")
library("tidyverse")

dat <- read.csv("Zayda_pcromp/alldataDavid.csv")

# Remove zeros
dat2 <- dat2[rowSums(dat2[,-1])!=0,]
# Make metadata file
meta <- data.frame(Sample.ID=colnames(dat2)[-1])
meta$Treatment <- gsub("_.*$","",meta$Sample.ID)
meta$Rep <- gsub("^.*_","",meta$Sample.ID)
meta$Colour <- rep(c("black","blue","darkred","purple","orange","darkgreen"), each=3)
meta$N2C3 <- meta$Treatment=="N2C3"

rownames(dat2) <- dat2$X
dat2_t <- t(as.matrix(dat2[,-1]))