data=readRDS('../../../data/GTEX/gtex_skin_count_data.rds')
col=readRDS('../../../data/GTEX/gtex_skin_colData.rds')
col2=col[col$gtex.smtsd=='Skin - Sun Exposed (Lower leg)',]
data2=data[,rownames(col2)]
colData=col2[col2$age_group!='middle',c('external_id','age_group')]
condition=factor(colData$age_group)
condition = relevel( condition, "young")
colData=data.frame(row.names=rownames(colData),condition)
countData=data2[,rownames(colData)]

###################################去除低表达和表达比例低的gene;并取整
library(edgeR)
countData <- countData[rowMeans(countData)>1,]
#keep <- filterByExpr(countData,min.count = 5)  
#countData=countData[keep,]
countData=round(countData)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
#dds$condition <- relevel(dds$condition, ref = "young")
########### 计算差异倍数
dds1 <- DESeq(dds) 

#将结果用result()函数来获取
#res <- lfcShrink(dds1, coef="condition_aged_vs_young", type="apeglm")
res <- results(dds1)
# res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
# 依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
res2=res1

# 获取padj（p值经过多重校验校正后的值）小于0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
res2_up<- res2[which(res2$log2FoldChange >0.58 & res2$padj < 0.05),]      # 表达量显著上升的基因
res2_down<- res2[which(res2$log2FoldChange <(-0.58) & res2$padj < 0.05),]    # 表达量显著下降的基因
res2_total <- rbind(res2_up,res2_down)

write.csv(res2_total,'dif.exposed.csv')
write.csv(res2_up,'up.exposed.csv')
write.csv(res2_down,'down.exposed.csv')
write.table(rownames(res2_down),'down.exposed.txt',row.names = F,quote = F,col.names = F)
write.table(rownames(res2_up),'up.exposed.txt',row.names = F,quote = F,col.names = F)

library(ggplot2)
library(ggrepel)
dat<-res2
dat$sig='none'
dat[rownames(res2_up),'sig']='up'
dat[rownames(res2_down),'sig']='down'
pdf("volcano_plot.exposed.pdf",height=12,width=11)
ggplot(dat,aes(x=log2FoldChange,y=-log10(padj),color=sig))+
geom_point()+
scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+  #确定点的颜色
theme_bw()+  #修改图片背景
theme(panel.border = element_rect(size = 2, fill = NA),axis.ticks.x=element_line(color="black",size=1.5),
      axis.ticks.y=element_line(color="black",size=1.5),axis.ticks.length.x = unit(0.3,'cm'), 
      axis.ticks.length.y = unit(0.3,'cm'),
  legend.title = element_blank()  #不显示图例标题
)+
theme(axis.title.x =element_text(size=14,face = "bold"), axis.title.y=element_text(size=14,face = "bold"),axis.text = element_text(size = 14,face = "bold")) +  #调整坐标轴字号
ylab('-log10 (p-adj)')+  #修改y轴名称
xlab('log2 (FoldChange)')+  #修改x轴名称
geom_vline(xintercept=c(-0.58,0.58),lty=3,col="black",lwd=0.5) +  #添加垂直阈值|FoldChange
geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.5)+  #添加水平阈值padj<0.05
#xlim(-max(dat$log2FoldChange)*1.2, max(dat$log2FoldChange)*1.2)
scale_x_continuous(limits=c(-max(dat$log2FoldChange)*1.2, max(dat$log2FoldChange)*1.2), breaks=seq(-3,3,1))
ggsave('volcano_plot.exposed.png',height=6,width=8)
dev.off()

# age rnk
colData=col2[,c('external_id','age2')]
condition=factor(colData$age2)
condition = relevel(condition, "20")
colData=data.frame(row.names=rownames(colData),condition)
countData=data2[,rownames(colData)]
###################################去除低表达和表达比例低的gene;并取整
library(edgeR)
countData <- countData[rowMeans(countData)>1,]
#keep <- filterByExpr(countData,min.count = 5)  
#countData=countData[keep,]
countData=round(countData)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
#dds$condition <- relevel(dds$condition, ref = "young")
########### 计算差异倍数
dds1 <- DESeq(dds) 
#将结果用result()函数来获取
#res <- lfcShrink(dds1, coef="condition_aged_vs_young", type="apeglm")
res=list()
for (i in 2:length(resultsNames(dds1))){
      resname=resultsNames(dds1)[i]
      tmp.res=results(dds1, name=resname)
    # res格式转化：用data.frame转化为表格形式
    res1 <- data.frame(tmp.res, stringsAsFactors = FALSE, check.names = FALSE)
    # 依次按照pvalue值log2FoldChange值进行排序
    res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
    res[[resname]]=res1
}
up=list()
down=list()
for(i in 1:length(res)){
    name=names(res)[i]
    print(name)
    tmp=res[[i]]
    rnk=data.frame(gene=rownames(tmp),lgfc=tmp$log2FoldChange)
    rnk=rnk[order(rnk$lgfc,decreasing = T),]
    #head(rnk,2)
    write.table(rnk,paste0(name,'.exposed.rnk'),row.names = F,quote = F,col.names = F,sep = '\t')
    tmp_up<- tmp[which(tmp$log2FoldChange >1 & tmp$padj < 0.05),] 
    tmp_down<- tmp[which(tmp$log2FoldChange <(-1) & tmp$padj < 0.05),] 
    up[[name]]=tmp_up
    down[[name]]=tmp_down
}



