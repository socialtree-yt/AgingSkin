rm(list = ls()) 
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
scenicLoomPath='sample_SCENIC.loom'
loom <- open_loom(scenicLoomPath)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)
saveRDS(regulons,'regulons.rds')
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)
sce=readRDS(file = "../05.2-subtype/scRNA_annodata.KC.rds")
table(Idents(sce))
sce=subset(sce,idents=c('SC','BC','SC_BC','GC_SC')) ###细胞数太少的cluster也不要
sce$celltype.g <- paste(Idents(sce), sce$group, sep = "_")
sce$celltype <- Idents(sce)
#Idents(sce) <- "celltype2"
Idents(sce) <- "celltype2"
head(sce@meta.data)

sub_regulonAUC <- regulonAUC[,match(colnames(sce),colnames(regulonAUC))]
dim(sub_regulonAUC)
sce 
#确认是否一致
identical(colnames(sub_regulonAUC), colnames(sce))

cellClusters <- data.frame(row.names = colnames(sce), 
                           seurat_clusters = as.character(sce$seurat_clusters))
cellTypes <- data.frame(row.names = colnames(sce), 
                        celltype = sce$celltype.g)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 
save(sub_regulonAUC,cellTypes,cellClusters,sce,file = 'for_rss_and_visual.Rdata')

