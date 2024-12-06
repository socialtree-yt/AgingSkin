######  step0 加载 各种R包  #####

rm(list=ls())
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
#library(scRNAseq)
if(!dir.exists('plot')){
    dir.create('plot')
}
load('for_rss_and_visual.Rdata')
head(cellTypes)
sub_regulonAUC[1:4,1:2] 
dim(sub_regulonAUC)
# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedResolution]) 

# 去除extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 
auc=getAUC(sub_regulonAUC) ########################################## 提取的regulon AUC 矩阵
dim(sub_regulonAUC)
#[1]  220 2638 #似乎没啥区别

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                    rowMeans(auc[,cells]))

# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
                                          center = T, scale=T)) 
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[]
regulonActivity_byGroup_Scaled=na.omit(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[,c('GC_SC_young','GC_SC_middle','GC_SC_aged',
                                                                 'SC_young','SC_middle','SC_aged',
                                                        'SC_BC_young','SC_BC_middle','SC_BC_aged','BC_young','BC_middle','BC_aged'
                                                        )]
head(regulonActivity_byGroup_Scaled,2)
library(pheatmap)
t=pheatmap(regulonActivity_byGroup_Scaled,cluster_cols = F)
g=rownames(regulonActivity_byGroup_Scaled)[t$tree_row$order]
ggsave('plot/complexheatmap.tree.scenic.png',t,width = 4.2,height = 9)
head(g,2)
rss <- calcRSS(AUC=auc, 
               cellAnnotation=cellTypes[colnames(sub_regulonAUC), 
                                           selectedResolution]) 
rss=rss[g,c('GC_SC_young','GC_SC_middle','GC_SC_aged',
                                                                 'SC_young','SC_middle','SC_aged',
                                                        'SC_BC_young','SC_BC_middle','SC_BC_aged','BC_young','BC_middle','BC_aged'
                                                        )]
head(rss,2)
library(tidyr)
to_long_data=function(cancer_count){
    cancer_count2=as.data.frame(t(as.matrix(cancer_count)))
    cancer_count2$sample=rownames(cancer_count2)
    cancer_count2=cancer_count2[,c('sample',rownames(cancer_count))]
    cancer_count2=cancer_count2%>% pivot_longer(.,-1,names_to = "gene",values_to = "expr") 
    return(cancer_count2)
}
data=as.data.frame(t(scale(t(regulonActivity_byGroup_Scaled[g,]))))
expr=as.data.frame(to_long_data(data))
expr$rss=''
for(i in 1:length(rownames(expr))){
    expr[i,'rss']=rss[expr[i,'gene'],expr[i,'sample']]
}
expr$rss=as.numeric(expr$rss)
expr$sample=factor(expr$sample,levels = c('GC_SC_young','GC_SC_middle','GC_SC_aged',
                                                                 'SC_young','SC_middle','SC_aged',
                                                        'SC_BC_young','SC_BC_middle','SC_BC_aged','BC_young','BC_middle','BC_aged'
                                                        ))
expr$gene=factor(expr$gene,levels = g)
ggplot(expr,aes(x=sample,y=gene))+
    geom_tile(mapping = aes(fill=expr))+
    geom_point(aes(size=rss))+#insert_left(phY, width = .09)+
    scale_fill_gradient2(low="#003366", high="#990033", mid="white")+scale_size(range=c(0, 2))+
        theme(panel.background = element_blank(),axis.line.x = element_line(color="black",size=0),
        axis.line.y = element_line(color="black",size=0),
        axis.text.x = element_text(size=5,angle = 90),axis.text.y = element_text(size=5),axis.line = element_line(colour = "black"),
             axis.ticks.x=element_line(color="black",size=0),legend.key = element_blank(),
      axis.ticks.y=element_line(color="black",size=0),axis.ticks.length.x = unit(0,'cm'), 
      axis.ticks.length.y = unit(0,'cm'))
ggsave('plot/complexheatmap.scenic.png',width = 4.2,height = 9)