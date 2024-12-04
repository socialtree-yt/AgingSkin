## heatmap
ECM=read.csv('ECM_matrisome_mouse.csv')
ECM=ECM[,c('Gene.Symbol','Division','Category','Synonyms','mouse')]
ECM=ECM[!duplicated(ECM$mouse),]
ECM=ECM[!is.na(ECM$mouse),]
rownames(ECM)=ECM$mouse
anno=ECM[bm.confirmed[bm.confirmed%in%ECM$mouse],]
#table(rownames(bm.confirmed)==bm.confirmed)

class=names(table(ECM$Category))

data=readRDS('../01DEG/fpkms.rds')
#data=data[,c('D30_Ctrl_6hr_TPA_Rep1','D30_Ctrl_6hr_TPA_Rep2','D30_IMQ_6hr_TPA_Rep1','D30_IMQ_6hr_TPA_Rep2')]
#rownames(data)=data$geneid

ex.col=data.frame(sample=colnames(data),group=c(rep('D6Ctrl',2),rep('D6IMQl',2),rep('D30Ctrl',2),
                                                                              rep('D30IMQ',2),rep('D30_Ctrl_6hr_TPA',2),rep('D30_IMQ_6hr_TPA',2)))
ex.col$age_group=ex.col$group

ex.sample=c('D30_Ctrl_6hr_TPA_Rep1','D30_Ctrl_6hr_TPA_Rep2','D30_IMQ_6hr_TPA_Rep1','D30_IMQ_6hr_TPA_Rep2')

all=unique(anno$mouse)
exp.all=data[all[all%in%rownames(data)],ex.sample]
exp.all=exp.all[rowMeans(exp.all)!=0,]
annocol=data.frame(row.names = colnames(exp.all),group=as.factor(c('D30_Ctrl_TPA','D30_Ctrl_TPA','D30_IMQ_TPA','D30_IMQ_TPA')))
annocolor=c("#00CCFF","#FF99CC")
names(annocolor)=c('D30_Ctrl_TPA','D30_IMQ_TPA')
plot2=pheatmap(exp.all,show_rownames = T,scale='row',clustering_method = 'ward.D', show_colnames = T,
               annotation_col = annocol,annotation_colors = list(group=annocolor),border_color = NA,
               cellwidth =40,cellheight = 7,
               #cellwidth =7,cellheight = 7,
                         cluster_cols = F,cluster_rows = T,#clustering_method ='centroid',
                         color=colorRampPalette(c("#0033CC","white","#FF3366"))(100), breaks=seq(-1.5, 1.5, length.out = 100)
    )
ggsave(paste0('BM','.fpkm.png'),plot2, width=12, height=12)


### customized GSEA
library(GseaVis)
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
getg=function(rnkpath){
    genes=read.table(rnkpath)
    gene=genes$V2
    names(gene)=genes$V1
    return(gene)
}
gene1=getg('CtlD30.vs.ImqD30.rnk')
gene2=getg('CtlD6.vs.ImqD6.rnk')
gene3=getg('CtlTPA.vs.ImqTPA.rnk')
pathway<-read.gmt("epiSC_age_vs_young.gmt")
y1 <- GSEA(gene1,TERM2GENE =pathway,pvalueCutoff = 1)
y2 <- GSEA(gene2,TERM2GENE =pathway,pvalueCutoff = 1)
y3 <- GSEA(gene3,TERM2GENE =pathway,pvalueCutoff = 1)

GSEAmultiGP(gsea_list = list(y1,y2,y3),geneSetID = "epiSC_age_vs_young_up",curve.col = ggsci::pal_lancet()(3),legend.position = "none",
            exp_name = c("CtlD30.vs.ImqD30","CtlD6.vs.ImqD6","CtlTPA.vs.ImqTPA"),addPval = T)+
            #exp_name = c("                   ","                 ","                     "),rect.bm.label = c('           ','              '))+
theme(axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_blank(),
        #axis.line.x=element_line(linetype=1,color="black",size=1),axis.line.y=element_line(linetype=1,color="black",size=1),
        #axis.ticks.x=element_line(color="black",size=1,lineend = 20),axis.ticks.y=element_line(color="black",size=1,lineend = 20),
        #axis.ticks.length = unit(0.4, "cm"),
        axis.text.x = element_blank(),axis.text.y = element_text(size=10),
        panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent'))
ggsave(paste0('epiSC_age_vs_young_up','.GSEA.png'), width=5, height=5)
