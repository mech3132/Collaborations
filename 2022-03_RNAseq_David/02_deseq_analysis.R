library(tidyverse)
library("ggVennDiagram")
library(network)
library(sna)
library(ggnet)


string_links <- read.delim("David_RNAseq/STRING/3702.protein.links.v11.5.txt", sep = " ")
string_links %>% separate(protein1, into=c(NULL,"protein1","b"), sep=".") %>%
  separate(protein2, into=c("a","protein2","b"), sep=".")
string_descrip <- read.delim("David_RNAseq/STRING/3702.protein.info.v11.5.txt")
