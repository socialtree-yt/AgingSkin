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
regulonActivity_byGroup_Scaled=regulonActivity_byGroup_Scaled[,c('GC_SC_aged','GC_SC_middle','GC_SC_young','SC_aged','SC_middle','SC_young',
                                                        'SC_BC_aged','SC_BC_middle','SC_BC_young','BC_aged','BC_middle','BC_young'
                                                        )]
                                                        ##rownames是矩阵列名，col是组名
col=data.frame(col=colnames(regulonActivity_byGroup_Scaled))
rownames(col)=col$col
sce=readRDS(file = "../05.2-subtype/scRNA_annodata.KC.rds")
table(Idents(sce))
ano_col=sce@meta.data[colnames(auc),c('group','celltype2')]
ano_col$celltype2=factor(ano_col$celltype2,levels = c('GC_SC','SC',
                                                        'SC_BC','BC'
                                                        ))
ano_col=ano_col[order(ano_col$celltype2),]
auc=auc[,rownames(ano_col)]
head(ano_col,2)
bc.col=ano_col[ano_col$celltype2=='BC',]
auc.bc=auc[,rownames(bc.col)]
y=auc.bc[,rownames(bc.col[bc.col$group=='young',])]
m=auc.bc[,rownames(bc.col[bc.col$group=='middle',])]
a=auc.bc[,rownames(bc.col[bc.col$group=='aged',])]
test=data.frame('gene'=rownames(y))
for (i in 1:length(rownames(y))){
    a1=a[i,]
    y1=y[i,]
    t=wilcox.test(as.numeric(a1),as.numeric(y1),alternative="two.sided",paired=F,var.equal=F,conf.level=0.95,exact=FALSE)
    #t=t.test(d1,d5,alternative="two.sided",paired=F,var.equal=F,conf.level=0.95)
    mean.a1=mean(a1)
    mean.y1=mean(y1)
    logfc=log2(mean.a1/mean.y1)
    test[i,'p_value']=t$p.value
    test[i,'logfc']=logfc
}
test=test[order(test$p_value),]
test$p.adj=p.adjust(test$p_value,method="bonferroni",n=length(test$p_value))
head(test,2)
res1=test
res1$p.adj[-log10(res1$p.adj)>100]=10^(-100)
# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res1_up<- res1[which(res1$logfc >0 & res1$p.adj < 0.05),]      # 表达量显著上升的基因
res1_down<- res1[which(res1$logfc <0 & res1$p.adj < 0.05),]    # 表达量显著下降的基因
res1_total <- rbind(res1_up,res1_down)
dim(res1_up)
dim(res1_down)
head(res1_total,2)
library(ggplot2)
library(ggrepel)
dat<-res1
dat$sig='none'
dat[rownames(res1_up),'sig']='up'
dat[rownames(res1_down),'sig']='down'
dat$sig=factor(dat$sig,levels=c('up','none','down'))
dat=dat[dat$logfc>(-1),]
pdf("plot/volcano_plot.auc.pdf",height=12,width=11)
ggplot(dat,aes(x=logfc,y=-log10(p.adj),color=sig))+
geom_point()+
scale_color_manual(values=c("#CC0000","#BBBBBB","#2f5688"))+  #确定点的颜色
theme_bw()+  #修改图片背景
theme(
  legend.title = element_blank()  #不显示图例标题
)+xlim(c(-1,1))+ylim(c(0,100))+
theme(axis.title.x =element_text(size=14,face = "bold"), axis.title.y=element_text(size=14,face = "bold"),
      axis.text = element_text(size = 14,face = "bold")) +  #调整坐标轴字号
ylab('-log10 (p-adj)')+  #修改y轴名称
xlab('log2 (FoldChange)')+  #修改x轴名称
#geom_vline(xintercept=c(-0.25,0.25),lty=3,col="black",lwd=0.5) +  #添加垂直阈值|FoldChange|>2
geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+  #添加水平阈值padj<0.05
geom_label_repel(data = dat, max.overlaps = 40,
                     aes(x = logfc, y = -log10(p.adj), label = gene),
                     size = 1.5,color="black",
                     #box.padding = unit(0.5, "lines"),
                     #point.padding = unit(0.8, "lines"), 
                     #segment.color = "black",   #连线的颜色
                     #segment.size = 0.4,  #连线的粗细
                     #arrow = arrow(length=unit(0.01, "npc")), #标签、点之间连线的箭头
                     show.legend = FALSE)
dev.off()
regulons=readRDS('regulons.rds')
head(regulons,2)
cellTypes2=sce@meta.data[,c('group','celltype2')]
colnames(cellTypes2)[2]='celltype'
head(cellTypes2)
selectedResolution2 <- "celltype"
rss2 <- calcRSS(AUC=auc, 
               cellAnnotation=cellTypes2[colnames(sub_regulonAUC), 
                                           selectedResolution2]) 
library(tidyr)
to_long_data=function(cancer_count){
    cancer_count2=as.data.frame(t(as.matrix(cancer_count)))
    cancer_count2$sample=rownames(cancer_count2)
    cancer_count2=cancer_count2[,c('sample',rownames(cancer_count))]
    cancer_count2=cancer_count2%>% pivot_longer(.,-1,names_to = "gene",values_to = "expr") 
    return(cancer_count2)
}
data2=to_long_data(rss2)
head(data2,2)
for(i in names(table(data2$sample))){
    tmp=data2[data2$sample==i,]
    tmp$rank=rank(tmp$expr)
    data2[data2$sample==i,'rank']=tmp$rank
}
head(data2,2)
regulon2=data2
regulon2=regulon2[order(regulon2$rank,decreasing = T),]
regu=function(gene){
    tmp=list(regulons[[gene]])
    return(tmp)
}
regulon2$target.gene=sapply(regulon2$gene,regu)
#for (i in 1:nrow(regulon1)) regulon1[i,'target.gene'] <- paste(regulon1[[i,'..target.gene..']], collapse = ',')
regulon2$target.gene=as.character(regulon2$target.gene)
#regulon1=unlist(regulon1)
head(regulon2,2)
library(dplyr)
data2%>%group_by(sample)%>%top_n(n=5,wt=rank)->top5.2
ggplot(data2, aes(rank, expr, color = sample)) +
  geom_line() +
  #scale_colour_manual(limits = c('a1', 'a2', 'a3', 'a4', 'a5', 'a6'), values = c('orange', 'purple', 'green', 'blue', 'red', 'gray40')) +
  labs(x = 'rank', y = 'Rss', color = NULL) +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), 
        legend.key = element_rect(fill = 'transparent'))+
ggrepel::geom_text_repel(inherit.aes = F,max.overlaps = 20, data = top5.2, aes(x = rank, y = expr, label = gene, color=sample))
ggsave('plot/rss.rank.percelltype.png',height=12,width=11)