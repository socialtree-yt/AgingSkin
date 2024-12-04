data=readRDS('../../../data/GTEX/gtex_skin_count_data.rds')
col=readRDS('../../../data/GTEX/gtex_skin_colData.rds')
col=col[col$gtex.smtsd=='Skin - Sun Exposed (Lower leg)',]
data=data[,rownames(col2)]
ex.col=col[col$group=='exposed',]
ex.col$age_group=paste0(ex.col$gtex.age,'_',ex.col$group)
ex.sample=rownames(ex.col)

sabm=read.csv('BMpipeline-SApaper.csv')
bm.confirmed=sabm[sabm$Category=='Basement membrane'&sabm$supportance=='confirmed','Gene.symbol']
ECM=read.csv('matrisome_hs_masterlist.csv')
ECM=ECM[,c('Gene.Symbol','Division','Category','Synonyms')]
rownames(ECM)=ECM$Gene.Symbol
anno=ECM[bm.confirmed[bm.confirmed%in%ECM$Gene.Symbol],]
#table(rownames(bm.confirmed)==bm.confirmed)

class=names(table(ECM$Category))
for(i in class){
    gene=ECM[ECM$Category==i,'Gene.Symbol']
    write.table(gene,paste0(i,'.txt'),row.names = F,quote = F,col.names = F)
}
write.table(unique(ECM$Gene.Symbol),paste0('allECM','.txt'),row.names = F,quote = F,col.names = F)

#### heatmap
for(i in class){
    gene=anno[anno$Category==i,'Gene.Symbol']
    plot.data=data[gene[gene%in%rownames(data)],ex.sample]
    
    cellTypes=ex.col
    selectedTypes='age_group'
    cellsPerGroup <- split(rownames(cellTypes), 
                       cellTypes[,selectedTypes]) 
    expr_byGroup <- sapply(cellsPerGroup,
                                  function(cells) 
                                     # apply(plot.data[,cells],1,median)
                                      rowMeans(plot.data[,cells])
                                    )
    expr_byGroup=as.data.frame(expr_byGroup)
    expr_byGroup=expr_byGroup[,c('20-29_exposed','30-39_exposed','40-49_exposed','50-59_exposed','60-69_exposed','70-79_exposed')]
    library(pheatmap)
    plot2=pheatmap(expr_byGroup,show_rownames = T,scale='row',clustering_method = 'ward.D', show_colnames = T,
               cellwidth =7,cellheight = 7,
                         cluster_cols = F,cluster_rows = T,#clustering_method ='centroid',
                         color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100)
    )
    ggsave(paste0(i,'.exposed.png'),plot2, width=3, height=8)
}


######  violin plot
to_long_data=function(cancer_count,groups){
    cancer_count2=as.data.frame(t(as.matrix(cancer_count)))
cancer_count2$sample=rownames(cancer_count2)
cancer_count2=cancer_count2[,c('sample',rownames(cancer_count))]
cancer_count2=cancer_count2%>% pivot_longer(.,-1,names_to = "gene",values_to = "expr") 
}
col.exposed=col[col$gtex.smtsd=='Skin - Sun Exposed (Lower leg)',]
data.exposed=data[,rownames(col.exposed)]
table(rownames(col.exposed)==colnames(data.exposed))
plot=function(data,col,dir){
    if(!dir.exists(dir)){
    dir.create(dir)
    }
    data2=to_long_data(data,'none')
    data2$age=col[data2$sample,'gtex.age']
    #head(data2)
    plot=data2[data2$gene%in%gene,]
   # p <- ggboxplot(plot, x = "gene", y = "expr",ylab = 'TPM',
    #          color = "age", palette = "nejm")+
    #theme(axis.title.x =element_text(size = 15),axis.title.y =element_text(size = 15))
    p=ggplot(data = plot,mapping = aes(x = gene,y = expr)) +
    geom_violin(scale = "width",adjust =1,trim = TRUE,mapping = aes(fill = gene)) +
    theme(panel.background = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),
       axis.text.x = element_text(size=10),axis.text.y = element_text(size=10),axis.line = element_line(colour = "black"))
    
    #theme(legend.position = 'top',axis.title = element_blank(),legend.text=element_text(margin = margin(r = 180, unit = "pt"),size=0),
    #      axis.text = element_blank(),
    #      axis.text.x = element_blank(), axis.text.y = element_blank(),legend.title=element_blank(),
    #      legend.key.size=unit(0.4,"inches"))+
    #guides(colour = guide_legend(override.aes = list(size=10)))
    #p + stat_compare_means(aes(group = group),label = "p.signif",size=20,hide.ns = TRUE,vjust =0.7)
    ggsave(filename = paste0(dir,'/plot.png'),plot = p,device = "png",width = 42,height = 15,units = "cm")
    #select=c('COL4A1','HSPG2','LAMC1','LAMA3','NID1','NID2')
    select=gene
    for(i in select){
        compar1=list(c('20-29','60-69'),c('20-29','70-79'))
        plot2=data2[data2$gene==i,]
        plot2$age=factor(plot2$age)
        Q <- quantile(plot2$expr, probs=c(.25, .75), na.rm = FALSE)
        iqr <- IQR(plot2$expr)
        up <- Q[2]+1.5*iqr # Upper Range 上限
    
        low<- Q[1]-1.5*iqr # Lower Range 下限
    
        #eliminated<- subset(plot2, plot2$expr > (Q[1] - 1.5*iqr) & plot2$expr < (Q[2]+1.5*iqr))
        p2=ggplot(data = plot2,mapping = aes(x = age,y = expr)) +labs(x = i,y='TPM')+ggtitle(dir)+
        geom_violin(scale = "width",adjust =1,trim = TRUE, fill = "steelblue")+
        geom_boxplot(width = 0.2)+#ylim(0,1)+  #记得删掉
        theme(panel.background = element_blank(),axis.line.x = element_line(color="black",size=0.8),axis.line.y = element_line(color="black",size=0.8),
        axis.text.x = element_text(size=10),axis.text.y = element_text(size=10),axis.line = element_line(colour = "black"),
             axis.ticks.x=element_line(color="black",size=0.8),
      axis.ticks.y=element_line(color="black",size=0.8),axis.ticks.length.x = unit(0.2,'cm'), 
      axis.ticks.length.y = unit(0.2,'cm'))+
        guides(colour = guide_legend(override.aes = list(size=10)))+
        stat_compare_means(comparisons = compar1,bracket.size = 0,method = 'wilcox.test',label = "p.signif",size=0,p.adjust.methods='bonferroni')
        #p2 <- ggboxplot(plot2, x = "age", y = "expr",ylab = '',xlab=i,bxp.errorbar = T)#+ coord_cartesian(ylim = c(0,up))
        ggsave(filename = paste0(dir,'/',i,'.png'),plot = p2,device = "png",width = 7,height = 10,units = "cm")
    }
}
gene=c('COL4A1','COL7A1','COL17A1','COL18A1')
plot(data.exposed,col.exposed,'tmp.ex')


#### dotplot
data=read.table('BM.GSEA.result.csv',sep=',',header = T)
library(ggplot2)
ggplot(data,aes(x = NES,y = group))+
    geom_point(aes(color = -log10(NOM.p.val+0.000001)),size=10)+
    #scale_color_gradient(low = "#FF9999", high = "#00CCFF", na.value = '#00CCFF',limits=c(0,0.05))+
    #scale_color_gradient(low = "blue", high = "red",breaks = c(1,-log10(0.05+0.000001),2,3,4,5,6))+
    scale_color_gradient(low = "blue", high = "red",breaks = c(1,-log10(0.05+0.000001),-log10(0.01+0.000001),-log10(0.001+0.000001),6))+
    #scale_colour_gradientn(colours=c( "orange","red","blue"),limits=c(0,0.1))+
    xlab("NES")+
    theme_bw(base_size = 20)+theme(panel.border = element_rect(size = 2, fill = NA),axis.ticks.x=element_line(color="black",size=1.5),
      axis.ticks.y=element_line(color="black",size=1.5),axis.ticks.length.x = unit(0.3,'cm'), 
      axis.ticks.length.y = unit(0.3,'cm'))+
    #edit legends
    guides(
        #reverse color order (higher value on top)
        color = guide_colorbar(reverse = TRUE))+
    xlim(-max(data$NES)*1.2, max(data$NES)*1.2)
ggsave('BM.GSEA.dotplot.png',width=10,height = 10)

