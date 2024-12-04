library(Seurat)
data=readRDS('../../02-cell_anotation/scRNA_annodata.rds')
sce=subset(data,idents = c('KC'))
DefaultAssay(sce) <- "RNA"
scRNAlist=SplitObject(sce,split.by = 'sample')
sce <- merge(scRNAlist[[1]], y=c(scRNAlist[[2]], scRNAlist[[3]], scRNAlist[[4]], scRNAlist[[5]], 
                                           scRNAlist[[6]], scRNAlist[[7]], scRNAlist[[8]], scRNAlist[[9]]))
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2500)
sce <- ScaleData(sce,verbose = FALSE)
sce <- RunPCA(sce,verbose = FALSE) 
sce=RunHarmony(sce, group.by.vars = "sample")
scRNA.integrated=sce
scRNA.integrated <- RunUMAP(object = scRNA.integrated, reduction = "harmony", dims = 1:19)
scRNA.integrated <- FindNeighbors(scRNA.integrated, reduction = "harmony", dims = 1:19, verbose=FALSE)
scRNA.integrated <- FindClusters(scRNA.integrated, reduction = "harmony", resolution = 2.8)  ##############################?????????????????
colors<-c("#4d648d", "#0450fb", "#11aac4", "#42e8f3", "#AEC7E8", "#2CA02C", "#98DF8A", "#9eccaf", "#daf400", "#983b59", "#e81f3f", "#ff8b94", "#ffd3b6", "#f9ae34", "#ffdb00", "#723584", "#9264eb", "#ff00ff", "#E377C2", "#de94e4", "#F7B6D2", "#C5B0D5", "#8C564B", "#C49C94", "#BCBD22", "#DBDB8D", "#7F7F7F", "#C7C7C7", "#a7a7b7", "#9999FF",  "#CC6633", "#990066", "#003333", "#996666","#CF6633")
DimPlot(scRNA.integrated, reduction = "umap", label = TRUE,group.by = 'group', repel = FALSE, cols=colors[1:length(unique(scRNA.integrated@active.ident))]) + NoLegend() + ggtitle("UMAP plot: colored by cluster")
DimPlot(scRNA.integrated, reduction = "umap", label = TRUE, cols=colors[1:length(unique(scRNA.integrated@active.ident))])
DimPlot(scRNA.integrated, reduction = "umap", label = TRUE,group.by = 'sample', repel = FALSE) + ggtitle("UMAP plot: colored by cluster")
saveRDS(scRNA.integrated, file = "subtype.harmony.rds")