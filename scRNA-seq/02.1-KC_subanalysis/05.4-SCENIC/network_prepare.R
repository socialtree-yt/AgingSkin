library(Seurat)
#library(SCopeLoomR)
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
#auc=readRDS(file = "auc.rds")
sce=readRDS(file = "../05.2-subtype/scRNA_annodata.KC.rds")
sce$celltype.g <- paste(sce$celltype2, sce$group, sep = "_")
head(sce,2)
scenicLoomPath='sample_SCENIC.loom'
loom <- open_loom(scenicLoomPath)

regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4] 
saveRDS(regulons_incidMat,'regulons_incidMat.rds')
load('for_rss_and_visual.Rdata')
head(cellTypes)
sub_regulonAUC[1:4,1:2] 
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
auc=getAUC(sub_regulonAUC)
head(auc,2)
bc_auc=auc[,rownames(cellTypes)[cellTypes$celltype%in%c('BC_young','BC_middle','BC_aged')]]
dim(bc_auc)
corr.result<-cor(t(bc_auc),method = 'pearson')
head(corr.result)
pheatmap(corr.result,scale='none')
ggsave('TF.interaction.cor.png',width = 13,height = 13)
#rownames(corr.result)[p2$tree_row$order]
regulons <- regulonsToGeneLists(regulons_incidMat)
head(regulons,1)
over=function(regulon,target){
    over=sum(target%in%regulon)
    return(over)
}
sabm=read.csv('BMpipeline-SApaper.csv')
head(sabm,2)
bm.confirmed=sabm[sabm$Category=='Basement membrane'&sabm$supportance=='confirmed','Gene.symbol']
length(bm.confirmed)
write.table(bm.confirmed,'bm.confirmed.txt',quote = F,row.names = F,col.names = F)
bm.pathway=unlist(readRDS('../05.3-pathscore/bm.pathway.rds')[c('basement membrane')])
length(bm.pathway)
tt2=sapply(regulons,over,target=unique(c(bm.confirmed,bm.pathway)))
tt2[order(tt2,decreasing = T)]
node=data.frame(TF=names(tt2),BMnumber=as.numeric(tt2))
head(node,2)
write.table(node,'node.txt',sep='\t',quote = F,row.names = F)
library(tidyr)
to_long_data=function(cancer_count){
    cancer_count2=as.data.frame(t(as.matrix(cancer_count)))
    cancer_count2$sample=rownames(cancer_count2)
    cancer_count2=cancer_count2[,c('sample',rownames(cancer_count))]
    cancer_count2=cancer_count2%>% pivot_longer(.,-1,names_to = "gene",values_to = "expr") 
    return(cancer_count2)
}
edge=to_long_data(corr.result)
colnames(edge)=c('from','to','weight')
edge$class='TF'
edge$fromsize=tt2[edge$from]
edge$tosize=tt2[edge$to]
edge=edge[edge$from!=edge$to,]
head(edge)
table(abs(edge$weight)>0.3)
edge2=edge[abs(edge$weight)>0.3,]
write.table(edge2,'edge.txt',sep='\t',quote = F,row.names = F)