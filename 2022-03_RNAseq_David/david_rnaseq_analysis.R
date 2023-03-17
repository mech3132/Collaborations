#!bin/bash

##### David analysis #######

dat <- read.delim("David_RNAseq/galaxy_data.txt")
str(dat)
dat$log2FC_abs <- abs(dat$log2.FC) # Get absolute value
dat$wald_abs = abs(dat$wald.statistic) # Get absolute value

## Wald ECDF

wald_ecdf <- ecdf(dat$wald_abs)
png("David_RNAseq/wald_cdf_all.png", height = 500, width = 700)
plot(sort(log10(dat$wald_abs)), wald_ecdf(sort(dat$wald_abs))
     , type = 'l'
     , xlab="Wald test statistic (abs log10)"
     , ylab = "Cumulative distirbution of values (cdf)"
     , main = "CDF of wald statistic, all data")
dev.off()

## Wald ECDF

FC_ecdf <- ecdf(dat$log2FC_abs)
png("David_RNAseq/FC_cdf_all.png", height = 500, width = 700)
plot(sort(log10(dat$log2FC_abs)), wald_ecdf(sort(dat$log2FC_abs))
     , type = 'l'
     , xlab="Wald test statistic (abs log10)"
     , ylab = "Cumulative distirbution of values (cdf)"
     , main = "CDF of wald statistic, all data")
dev.off()

# is FC strongly correlated with statistical significance?
png("David_RNAseq/corr_FC_wald.png")
plot(log(dat$log2FC_abs), log(dat$wald_abs)
     , xlab="Log2 FC (abs log)"
     , ylab="Wald statistic (abs log)")
dev.off()
# Yes

### Try looking at loadings

dat2 <- read.csv("Zayda_pcromp/alldataDavid.csv")
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

pcoa_dat2t <- prcomp(dat2_t)
var_explained <- pcoa_dat2t$sdev^2/sum(pcoa_dat2t$sdev^2)

png(filename = "David_RNAseq/pcoa_allsamples.png")
plot(pcoa_dat2t$x[,1], pcoa_dat2t$x[,2]
     , xlab=paste0("PC1 (",round(var_explained[1]*100,2),"% variation explained)")
     , ylab=paste0("PC2 (",round(var_explained[2]*100,2),"% variation explained)")
     , col = meta[match(rownames(pcoa_dat2t$x), meta$Sample.ID),"Colour"]
     , pch = 19)
legend("bottomright",legend=unique(meta$Treatment), col=unique(meta$Colour), pch=19, border = NULL, cex = 1)
dev.off()
# Get loading of each gene for PC1; and then see whether PC1 is the BEST option for each gene or not.
PC_rankings <- t(apply(abs(pcoa_dat2t$rotation), MARGIN = 1, function(x) x/max(x)))
PC_perctotalrank <- t(apply(abs(pcoa_dat2t$rotation), MARGIN = 1, function(x) x/sum(x)))

PC1_corr <- data.frame(gene = rownames(pcoa_dat2t$rotation), PC_loadings = pcoa_dat2t$rotation[,1], PC_rankings = PC_rankings[,1], PC_proprank = PC_perctotalrank[,1])
# PC1_corr_filt <- PC1_corr[-which(PC1_corr$PC_loadings==0),]

par(mar=c(5.1, 5.1, 4.1, 2.1))
plot((PC1_corr$PC_loadings), PC1_corr$PC_proprank
     , xlab="PC1 loadings for each gene"
     , ylab="Relative loading\ncompared to other PCs (%)"
     )

####### Heatmap
cor(dat2_t)



####### Try ANCOM on filtered dataset
source("../Code/ANCOM_updated_code_MYC.R")
# Remove things that are very low read count-- aka, things that are less than 1 normalized read
dat2_t_filt <- dat2_t[,colSums(dat2_t)>0]
# There are still an insane number of sequences... let's look at variance in dataset 
hist(log(apply(dat2_t_filt, MARGIN=2, function(x) sd(x))), breaks=1000)
# Filter by sequences whose mean-standardized SD is > 1
meanstvar <- apply(dat2_t_filt, MARGIN=2, function(x) sd(x)/mean(x))
genes_to_keep <- names(which(meanstvar>1))

# Convert to ANCOM format
ancom_data <- data.frame(Sample.ID = rownames(dat2_t_filt), exp(dat2_t_filt))

# Filter data
# tempgene_keep <- names(sort(-abs(colSums(ancom_data[,-1])))[1:100])
ancom_data_subset <- ancom_data[,c("Sample.ID", genes_to_keep)]

ancom_rnaseq <- ANCOM.main.myc(OTUdat=ancom_data_subset, Vardat = meta, main.var = "N2C3"
                       , repeated = FALSE # Repated measures; default false
                       , adjusted = FALSE # Additional factors
                       , multcorr = 1 # strigency of corrections; 1 is conservattive, 3 is no correction
                       , prev.cut = 0.9 # prevalence threshold
                       , sig = 0.05 # significance threshold for KW tests
                       )

ancom_rnaseq$Plot.volcano
ancom_rnaseq$W.taxa[ancom_rnaseq$W.taxa$detected_wcrit,]
