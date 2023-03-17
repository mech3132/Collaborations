#!/bin/bash

conda activate qiime2-2022.8

# Remove labels in ESVs

cat repset_16S_tim.fasta | sed 's/ Bacteria.*//g' | sed 's/ Archaea.*//g' | sed 's/ Eukaryota.*//g' | sed 's/ NA.*//g' > repset_16S_tim_cleaned.fasta

qiime tools import \
--input-path repset_16S_tim_cleaned.fasta \
--output-path repset_16S_tim_cleaned.qza \
--type FeatureData[Sequence]

# Clean up
rm repset_16S_tim_cleaned.fasta

## Run the tree insertion on compute canada
# Requires ~ 22h with 16G RAM per node on 10 nodes
wget 'https://data.qiime2.org/2022.8/common/sepp-refs-silva-128.qza'
mkdir -p tim_sepp_tree_16S
 qiime fragment-insertion sepp \
  --i-representative-sequences repset_16S_tim_cleaned.qza \
  --i-reference-database sepp-refs-silva-128.qza \
  --p-threads 10 \
  --o-tree tim_sepp_tree_16S/insertion-tree_16S_tim.qza \
  --o-placements tim_sepp_tree_16S/insertion-placements_16S_tim.qza
  
 qiime tools export \
 --input-path tim_sepp_tree_16S/insertion-tree_16S_tim.qza \
 --output-path tim_sepp_tree_16S/insertion-tree_16S_silva128
 
  qiime tools export \
 --input-path tim_sepp_tree_16S/insertion-placements_16S_tim.qza \
 --output-path tim_sepp_tree_16S/insertion-tree_16S_silva128
 
 # Make list of ESVs to keep
 echo -e "OTUID" > mock_otu.txt
 grep ">" repset_16S_tim_cleaned.fasta | sed 's/>//g' >> mock_otu.txt
 
 # Run R code
 # First, need to expand size of membory R is allowed to use
 ulimit -s 16384
 # Then you can run in R
 Rscript filter_Tree.R
 
 # Move trees to easily accesible place
 mkdir -p SEPP_TREE
 mv tim_sepp_tree_16S/insertion-tree_16S_silva128/tree.nwk SEPP_TREE/
 mv tim_sepp_tree_16S/insertion-tree_16S_silva128/filtered-tree.nwk SEPP_TREE/
 
 

