library(ape)
library(tidyverse)

tre <- read.tree("tim_sepp_tree_16S/insertion-tree_16S_silva128/tree.nwk")
dat <- read.delim("mock_otu.txt")
dat_filt <- dat %>% filter(OTUID %in% tre$tip.label) 
tre_filt <- keep.tip(tre, dat_filt$OTUID)
write.tree(tre_filt, file="tim_sepp_tree_16S/insertion-tree_16S_silva128/filtered-tree.nwk")
