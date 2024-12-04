library(Seurat)
library(ggpubr)
library(ggsignif)
library(ggplot2)
obj.integrated=readRDS(file = "scRNA_annodata.KC.BMscore.rds")
#obj.integrated=subset(obj.integrated,idents = c('SC','BC','SC_BC'))
Idents(obj.integrated) <- "celltype2" #######也就是合并的celltype
obj.integrated=subset(obj.integrated,idents = c('SC','SC_BC','BC'))

obj.integrated$celltype.g <- paste(Idents(obj.integrated), obj.integrated$group, sep = "_")
Idents(obj.integrated) <- "celltype.g"
gene=readRDS('bm.pathway.rds')
exp=t(as.matrix(obj.integrated@assays$RNA@data))
if(!dir.exists('corplot-BM')){
    dir.create('corplot-BM')
}
library(psych)
g=unique(gene$`basement membrane`)
g=g[g%in%colnames(exp)]
cor_matrix=data.frame(gene=g,r=0,p=0)
rownames(cor_matrix)=cor_matrix$gene
myCor = corr.test(exp[,'COL17A1'],meta[,'basement.membrane1'],  
               use = "pairwise", # 缺失值处理的方式
               method="pearson", # 计算相关性的方法有"pearson", "spearman", "kendall"
               adjust = "none"   # p值矫正的方法
)
for(i in g){
    tt=as.data.frame(cbind(exp[,i],meta[,'basement.membrane1']))
    colnames(tt)=c(i,'basement.membrane1')
    myCor = corr.test(tt[,i],tt[,'basement.membrane1'],  
               use = "pairwise", # 缺失值处理的方式
               method="pearson", # 计算相关性的方法有"pearson", "spearman", "kendall"
               adjust = "none"   # p值矫正的方法
    )
    cor_matrix[i,'r']=myCor$r
    cor_matrix[i,'p']=myCor$p
    plt=ggscatter(tt, x = i, y = "basement.membrane1",color = "#00AFBB",size=1,
          add = "reg.line", conf.int = TRUE,    
          add.params = list(color = "red", fill = "lightgray")
          
    )+
    stat_cor(method = "pearson")
    ggsave(paste0('corplot-BM/',i,'.png'),plt,width = 5,height = 5)
}
#cor_matrix=cor_matrix[cor_matrix$p<0.05,]
cor_matrix=cor_matrix[order(cor_matrix$r,decreasing = T),]
cor_matrix$rank=length(cor_matrix$gene):1
head(cor_matrix,10)
top10=cor_matrix[1:10,]
top10$gene
ggplot(cor_matrix, aes(x=rank, y=r))+
  geom_point(color='#9900CC',size=4)+
  #geom_hline(yintercept = c(2,-2), linetype=2, size=0.25)+
  #geom_hline(yintercept = c(0), linetype=1, size=0.5)+
  #geom_vline(xintercept = sum(pre_ranked_all_genes$trend == 'down')+0.5, linetype=2, size=0.25)+
  #ggrepel::geom_text_repel(inherit.aes = F, data = top10, aes(x=rank, y=r, label=gene, color='red'), 
  #                         size=3, direction = 'y')+
  #scale_color_manual(values=c('bottom5'='blue', 'top5'='red', 'medium'='black'))+
  scale_size_continuous(name='-log10(FDR)')+
  #scale_y_continuous(breaks=c(-5, -2, 0, 2, 5, 10))+
  xlab('rank of differentially expressed genes')+
  #theme_bw()+
  theme(panel.background = element_blank(),axis.line.x = element_line(color="black",size=1.8),
        axis.line.y = element_line(color="black",size=1.8),
        axis.text.x = element_text(size=0),axis.text.y = element_text(size=0),axis.line = element_line(colour = "black"),
             axis.ticks.x=element_line(color="black",size=1.8),legend.key = element_blank(),
      axis.ticks.y=element_line(color="black",size=1.8),axis.ticks.length.x = unit(0.5,'cm'), 
      axis.ticks.length.y = unit(0.5,'cm'))
ggsave('corplot-BM/BMrankplot.png',width = 11,height = 10)
