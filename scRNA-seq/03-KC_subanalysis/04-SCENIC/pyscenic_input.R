### Initialize settings
library(Seurat) 
data=readRDS(file = "../05.2-subtype/scRNA_annodata.KC.rds")
data=subset(data,idents=c('SC','BC','SC_BC','GC_SC')) ###细胞数太少的cluster也不要
t=t(as.matrix(data@assays$RNA@counts))
head(t)
dim(t)
write.csv(t(as.matrix(data@assays$RNA@counts)),file = "RNA_count.csv")