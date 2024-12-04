library(Seurat)
library(ggpubr)
library(ggsignif)
library(ggplot2)
obj.integrated=readRDS(file = "../05.2-subtype/scRNA_annodata.KC.rds")
#obj.integrated=subset(obj.integrated,idents = c('SC','BC','SC_BC'))
Idents(obj.integrated) <- "celltype2" #######也就是合并的celltype
obj.integrated=subset(obj.integrated,idents = c('GC_SC','SC','SC_BC','BC'))
obj.integrated2=subset(obj.integrated,idents = c('BC'))
obj.integrated$celltype.g <- paste(Idents(obj.integrated), obj.integrated$group, sep = "_")
Idents(obj.integrated) <- "celltype.g"
table(Idents(obj.integrated))
gene=readRDS('bm.pathway.rds')
obj.integrated=AddModuleScore(obj.integrated,features = list(gene[[1]]),name=names(gene)[1])
options(repr.plot.width=12, repr.plot.height=12)
Idents(obj.integrated)=obj.integrated$celltype.g
DefaultAssay(obj.integrated)='RNA'
levels(obj.integrated)=c('GC_SC_young','GC_SC_middle','GC_SC_aged','SC_aged','SC_middle','SC_young','SC_BC_aged',
                         'SC_BC_middle','SC_BC_young','BC_aged','BC_middle','BC_young')
VlnPlot(obj.integrated,features = 'basement.membrane1')
head(obj.integrated,2)
gene1=c('basement.membrane1')
exprs <- data.frame(FetchData(object = obj.integrated, vars = c(gene1,'celltype2','celltype.g')))
saveRDS(obj.integrated,'scRNA_annodata.KC.BMscore.rds')
noise <- rnorm(n = length(x = exprs[,gene1])) / 100000
exprs[,gene1] <- exprs[, gene1] + noise
#exprs$celltype.g=factor(exprs$celltype.g,levels = c('SC_aged','SC_middle','SC_young','SC_BC_aged','SC_BC_middle','SC_BC_young','BC_aged','BC_middle','BC_young'))
exprs$celltype.g=factor(exprs$celltype.g,levels = c('GC_SC_young','GC_SC_middle','GC_SC_aged','SC_young','SC_middle','SC_aged',
                                                    'SC_BC_young','SC_BC_middle','SC_BC_aged',
                                                  'BC_young','BC_middle','BC_aged'))
#gene1=c('basement.membrane1','basement.membrane.disassembly1','basement.membrane.organization1','basement.membrane.assembly1',
#       'basement.membrane.collagen.trimer1','regulation.of.basement.membrane.organization1')
if(!dir.exists('aging')){
    dir.create('aging')
}
for(i in gene1){
    #compar=list(c('BC_aged','BC_middle'),c('BC_middle','BC_young'),c('BC_aged','BC_young'))
    compar=list(c('BC_aged','BC_middle'),c('BC_aged','BC_young'),c('BC_middle','BC_young'))
    tmp=exprs[,c('celltype.g',i)] ########取出这个通路score，防止ggplot里y输入是字符导致无法计算
    colnames(tmp)=c('celltype.g','pathway')
    plot=ggplot(data = tmp,mapping = aes(x = celltype.g,y = pathway)) +
    geom_violin(scale = "width",adjust =1,trim = TRUE,mapping = aes(fill = celltype.g)) +
    geom_boxplot(width=.2,col="black",fill="white")+
    theme(panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
       axis.text.x = element_text(size=10),axis.text.y = element_text(size=10),axis.line = element_line(colour = "black"))+ggtitle(i)+
    stat_compare_means(comparisons = compar,label = "p.signif",method = 't.test',size=10,p.adjust.methods='bonferroni') #,method = 't.test'
    ggsave(paste0('aging/',i,'.pdf'),plot,width = 10,height = 10)
    ggsave(paste0('aging/',i,'.png'),plot,width = 10,height = 10)
}

noise <- rnorm(n = length(x = exprs[,gene1])) / 100000
exprs[,gene1] <- exprs[, gene1] + noise
#exprs$celltype.g=factor(exprs$celltype.g,levels = c('SC_aged','SC_middle','SC_young','SC_BC_aged','SC_BC_middle','SC_BC_young','BC_aged','BC_middle','BC_young'))
exprs$celltype.g=factor(exprs$celltype.g,levels = c('GC_SC_young','GC_SC_middle','GC_SC_aged','SC_young','SC_middle','SC_aged',
                                                    'SC_BC_young','SC_BC_middle','SC_BC_aged',
                                                  'BC_young','BC_middle','BC_aged'))
#gene1=c('basement.membrane1','basement.membrane.disassembly1','basement.membrane.organization1','basement.membrane.assembly1',
#       'basement.membrane.collagen.trimer1','regulation.of.basement.membrane.organization1')
if(!dir.exists('aging')){
    dir.create('aging')
}
for(i in gene1){
    #compar=list(c('BC_aged','BC_middle'),c('BC_middle','BC_young'),c('BC_aged','BC_young'))
    compar=list(c('BC_aged','BC_middle'),c('BC_middle','BC_young'),c('BC_aged','BC_young'))
    tmp=exprs[,c('celltype.g',i)] ########取出这个通路score，防止ggplot里y输入是字符导致无法计算
    colnames(tmp)=c('celltype.g','pathway')
    plot=ggplot(data = tmp,mapping = aes(x = celltype.g,y = pathway)) +
    geom_violin(scale = "width",adjust =1,trim = TRUE,mapping = aes(fill = celltype.g)) +
    geom_boxplot(width=.2,col="black",fill="white")+
    theme(panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
          axis.line.x=element_line(linetype=1,color="black",size=1),axis.line.y=element_line(linetype=1,color="black",size=1),
        axis.ticks.x=element_line(color="black",size=1,lineend = 20),axis.ticks.y=element_line(color="black",size=1,lineend = 20),
        axis.ticks.length.y = unit(0.4, "cm"),axis.ticks.length.x = unit(0, "cm"),
       axis.text.x = element_text(size=0),axis.text.y = element_text(size=0),axis.line = element_line(colour = "black"))+#+ggtitle(i)+
    stat_compare_means(comparisons = compar,label = "p.signif",bracket.size = 0,size=0,p.adjust.methods='bonferroni')
    #ggsave(paste0(i,'.pdf'),plot,width = 10,height = 10)
    ggsave(paste0('aging/',i,'.figure.png'),plot,width = 10,height = 10)
}
