#!bin/bash

all <- read.delim("alldataDavid.csv", sep=",")
# Made rownames == X
rownames(all) <- all[,1]
# transpose everything except X
allt <- t(all[,-1])

## Resume code pca ZMM
res.pca <- prcomp(allt)

fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     
)
ind <- get_pca_ind(res.pca)
ind
head(ind$coord)
v <- data.frame(ind$coord, max.print = 200)
v
write.csv(v,file="coordinatesflg22_alldavid.csv",row.names=F)