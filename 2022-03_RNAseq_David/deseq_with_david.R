# For DESEQ only

library("tidyverse")
library("DESeq2")
dat <- read.delim("David_RNAseq/raw_counts.txt")

# Remove gene lengths
dat <- dat[,-c(which(colnames(dat)=="gene_lengths_bp"))]
# Remove zeros
dat <- dat[rowSums(dat[,-1])!=0,]
rownames(dat) <- dat$Geneid
dat <- dat[,-1]

# Make metadata file
meta <- data.frame(SampleID=colnames(dat)) %>%
  separate(SampleID, sep = "_", into=c("Treatment","Rep"))
rownames(meta) <- colnames(dat)
# Make into deseq object
dds <- DESeqDataSetFromMatrix(countData = dat,
                              colData = meta,
                              design = ~ Treatment)
# Make MgSO4 first; control
dds$Treatment <- relevel(dds$Treatment, "MgSO4")

deseqRun <- DESeq(dds)
# 
# N2C3_res <- data.frame(results(deseqRun, contrast=c("Treatment","N2C3", "MgSO4")), Treatment="N2C3")%>%
#   rownames_to_column(var="GeneID")
# WCS365_res <- data.frame(results(deseqRun, contrast=c("Treatment","WCS365", "MgSO4")), Treatment="WCS365")%>%
#   rownames_to_column(var="GeneID")
# flg22_res <- data.frame(results(deseqRun, contrast=c("Treatment","flg22", "MgSO4")), Treatment="flg22")%>%
#   rownames_to_column(var="GeneID")
# pep1_res <- data.frame(results(deseqRun, contrast=c("Treatment","pep1", "MgSO4")), Treatment="pep1")%>%
#   rownames_to_column(var="GeneID")
# syrsyp_res <- data.frame(results(deseqRun, contrast=c("Treatment","syrsyp", "MgSO4")), Treatment="syrsyp")%>%
#   rownames_to_column(var="GeneID")
allRes <- rbind(data.frame(results(deseqRun, contrast=c("Treatment","N2C3", "MgSO4")), Treatment="N2C3")%>%
                  rownames_to_column(var="GeneID")
                ,data.frame(results(deseqRun, contrast=c("Treatment","WCS365", "MgSO4")), Treatment="WCS365")%>%
                  rownames_to_column(var="GeneID")
                , data.frame(results(deseqRun, contrast=c("Treatment","flg22", "MgSO4")), Treatment="flg22")%>%
                  rownames_to_column(var="GeneID")
                ,data.frame(results(deseqRun, contrast=c("Treatment","pep1", "MgSO4")), Treatment="pep1")%>%
                  rownames_to_column(var="GeneID")
                ,data.frame(results(deseqRun, contrast=c("Treatment","syrsyp", "MgSO4")), Treatment="syrsyp")%>%
                  rownames_to_column(var="GeneID")) 

# Total expression of genes in dataset
hist(log10(unlist(dat)), breaks=100) # ALL reads
# Summarize by gene mean
hist(log10(rowMeans(dat)), breaks=100)
abline(v=1.5, col="red", lty=2)

# Remove genes whose row means have expression < 10 counts
low_expression_genes <- names(which(log10(rowMeans(dat))<1.5))

# Filter Res genes by excluding low expression genes
allRes_filt <- allRes %>% filter(!(GeneID %in% low_expression_genes))

# Make table of just padj
allRes_filt_padj <- allRes_filt %>% 
  select(GeneID, Treatment, padj) %>%
  pivot_wider( names_from=Treatment, values_from=padj ) %>% as.data.frame()
rownames(allRes_filt_padj) <- allRes_filt_padj$GeneID
allRes_filt_padj <- allRes_filt_padj[,-which(colnames(allRes_filt_padj)=="GeneID")]

## Probe different p-value cutoffs to see how the proportion of shared genes
# change when you use different thresholds
# Filter dataset by those that are GREATER than threshold (not significant), and calculate number
# of treatments it is also NOT significant in. Higher threshold (0.05>0.01) should result in less shared genes
pval_cutoff=10^seq(log10(10^(-5)),log10(1), length.out=100)
sharedGenes <- t(apply(data.frame(cutoff=pval_cutoff), MARGIN=1
              , function(x) c(mean=mean(rowSums(allRes_filt_padj>x), na.rm=TRUE)
                              ,sd=sd(rowSums(allRes_filt_padj>x), na.rm=TRUE))
              ))
sharedGenesSummary <- data.frame(pval_cutoff,sharedGenes)
ggplot(sharedGenesSummary) + 
  geom_line(aes(x=pval_cutoff, y=mean)) +
  geom_errorbar(aes(x=pval_cutoff, ymin=mean-sd, ymax=mean+sd)) +
  scale_x_log10()+
  geom_vline(aes(xintercept=0.05), col="red")
# A threshold of 0.05 actually looks pretty reasonable

##### VENN DIAGRAM #######

# The number of mock-common-genes shared across treatments-- "core" should be the very centre
sharedwithmock <- list(
  N2C3=allRes_filt %>%
    filter(padj>0.05, Treatment=="N2C3") %>% pull(GeneID)
  ,WCS365=allRes_filt %>%
    filter(padj>0.05, Treatment=="WCS365") %>% pull(GeneID)
  ,flg22=allRes_filt %>%
    filter(padj>0.05, Treatment=="flg22") %>% pull(GeneID)
  ,syrsyp=allRes_filt %>%
    filter(padj>0.05, Treatment=="syrsyp") %>% pull(GeneID)
  ,pep1=allRes_filt %>%
    filter(padj>0.05, Treatment=="pep1") %>% pull(GeneID)
)
ggVennDiagram(sharedwithmock)

## Get core
core_across_all <- Reduce(intersect,sharedwithmock[c("N2C3","WCS365","flg22","pep1","syrsyp")] )

# What are the 33% that are in everything EXCEPT N2C3??
pooled_nonN2C3_sharedwithmock <- Reduce(intersect,sharedwithmock[c("WCS365","flg22","pep1","syrsyp")] )
N2C3_diff_from_core <- pooled_nonN2C3_sharedwithmock[which(!pooled_nonN2C3_sharedwithmock%in% sharedwithmock$N2C3)]
## list of genes that are common between mock and all other treatments, but not N2C3
length(N2C3_diff_from_core)

## Filter by genes that are DIFFERENT from mock-- how does differential expression compare across treatments?
difffrommock <- list(
  N2C3=allRes_filt %>%
    filter(padj<0.05, Treatment=="N2C3") %>% pull(GeneID)
  ,WCS365=allRes_filt %>%
    filter(padj<0.05, Treatment=="WCS365") %>% pull(GeneID)
  ,flg22=allRes_filt %>%
    filter(padj<0.05, Treatment=="flg22") %>% pull(GeneID)
  ,syrsyp=allRes_filt %>%
    filter(padj<0.05, Treatment=="syrsyp") %>% pull(GeneID)
  ,pep1=allRes_filt %>%
    filter(padj<0.05, Treatment=="pep1") %>% pull(GeneID)
)

ggVennDiagram(difffrommock)

# Get the 6133 genes that are unique to N2C3-- same as the other ones??
all_difffrommock_exceptN2C3 <- unique(unlist(difffrommock[c("WCS365","flg22","pep1","syrsyp")]))
diffmock_uniquetoN2C3 <- difffrommock$N2C3[which(!difffrommock$N2C3 %in% all_difffrommock_exceptN2C3)]

# Make network
diffmock_uniquetoN2C3
string_links

# Get genes that are common between different combos of N2C3/pep1/syrsyp
allRes_filt_padj
allRes_filt %>% filter(padj<0.05, Treatment %in% c("N2C3","pep1","syrsyp"))

# Loop through all treatments except MgSO4
for ( i in unique(meta$Treatment)[unique(meta$Treatment)!="MgSO4"] ) {
  tempMeta <- meta %>% filter(Treatment %in% c("MgSO4",i)) 
  tempDat <- dat %>% dplyr::select(one_of(tempMeta$SampleID))
  
  #dirtytransform
  dds <- DESeqDataSetFromMatrix(countData = tempDat,
                                colData = tempMeta,
                                design = ~ Treatment)
  deseqrun <- DESeq(dds)
  res <- results(deseqrun)
  res_tab <- data.frame(compare=i
                        , baseMean=res$baseMean
                        , log2FoldChange=res$log2FoldChange
                        , lfcSE=res$lfcSE
                        , stat=res$stat
                        , pvalue=res$pvalue
                        , padj=res$padj)
  rownames(res_tab) <- rownames(N2C3_dat_deseq)
}



###### scratch

N2C3_meta <- meta %>% filter(Treatment %in% c("MgSO4","N2C3")) 
N2C3_dat <- dat %>% dplyr::select(one_of(N2C3_meta$SampleID))

#dirtytransform
dds <- DESeqDataSetFromMatrix(countData = N2C3_meta,
                              colData = N2C3_meta,
                              design = ~ Treatment)
dds
deseqrun <- DESeq(dds)
res <- results(deseqrun)

res_tab <- data.frame(baseMean=res$baseMean
      , log2FoldChange=res$log2FoldChange
      , lfcSE=res$lfcSE
      , stat=res$stat
      , pvalue=res$pvalue
      , padj=res$padj)
rownames(res_tab) <- rownames(N2C3_dat_deseq)

ggplot(res_tab) +
  geom_point(aes(x=log2FoldChange, y=abs(stat), col=pvalue>0.05))


## Find all apirwise combos; shee shared non different.