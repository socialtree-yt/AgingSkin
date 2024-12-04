getfpkm=function(sample.dir){
    samples=read.table(paste0(sample.dir,'/samples.lst'),fill = T)
file=samples$V1
fpkm=read.table(paste0(sample.dir,'/',file[1],'.count'),header = T)[,c(1,7)]
colnames(fpkm)[2]=file[1]
fpkm=fpkm[order(fpkm[,1]),]
for(i in file[2:length(file)]){
    #if(i=='08_deep_senescence_rep2_SRR_029'){i='08_deep_senescence_rep2'}
    tt=read.table(paste0(sample.dir,'/',i,'.count'),header = T)[,c(1,7)]
    colnames(tt)[2]=i
    tt=tt[order(tt[,1]),]
    #if(i=='08_deep_senescence_rep2'){colnames(tt)[2]='08_deep_senescence_rep2_SRR_029'}
    fpkm=cbind(fpkm,tt)
    #fpkm=merge(fpkm,tt,by = 'gene_short_name')
}
fpkm=fpkm[,c('Geneid',file)]
    return(fpkm)
}
fpkm=getfpkm('AgingSkinTissue/')

library(edgeR)
fpkm=avereps(fpkm,fpkm$Geneid)
#head(fpkm)
rownames(fpkm)=fpkm[,1]
fpkm=fpkm[,-1]
fpkm=as.data.frame(fpkm)
name=rownames(fpkm)
fpkm=apply(fpkm,2,as.numeric)
rownames(fpkm)=name

saveRDS(fpkm,'counts.rds')


## calculate degs
fpkm=readRDS('counts.rds')

col=data.frame(row.names = colnames(fpkm),'external_id'=colnames(fpkm),age2=c(rep('old',3),rep('young',3)))
#head(col,2)
countData=fpkm
colData=col
condition=factor(colData$age2)
condition = relevel(condition, "Ctl")
colData=data.frame(row.names=rownames(colData),condition)

table(rownames(colData)==colnames(countData))
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
    tmp.res=results(dds1,pAdjustMethod = 'BH', name=resname)
    # res格式转化：用data.frame转化为表格形式
    res1 <- data.frame(tmp.res, stringsAsFactors = FALSE, check.names = FALSE)
    # 依次按照pvalue值log2FoldChange值进行排序
    res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
    res[[resname]]=res1
}
normfpkm=counts(dds1, normalized=TRUE)
normfpkm=as.data.frame(normfpkm)

write.table(normfpkm,'norm.exp.count.txt',row.names = T,quote = F,col.names = T,sep='\t')
up=list()
down=list()
for(i in 1:length(res)){
    name=names(res)[i]
    print(name)
    tmp=res[[i]]
    tmp_up<- tmp[which(tmp$log2FoldChange >1 ),] # & tmp$padj < 0.05
    tmp_down<- tmp[which(tmp$log2FoldChange <(-1)),] 
    up[[name]]=tmp_up
    down[[name]]=tmp_down
}
up.oy.gene=unique(rownames(up[['condition_old_vs_young']]))
down.oy.gene=unique(rownames(down[['condition_old_vs_young']]))
write.table(up.oy.gene,'up.oy.gene.count.txt',row.names = F,quote = F,col.names = F)
write.table(down.oy.gene,'down.oy.gene.count.txt',row.names = F,quote = F,col.names = F)


## get rnk
rnkm=function(deg){
    tt=data.frame('gene'=rownames(deg),'logfc'=deg$log2FoldChange)
    tt=tt[order(tt$logfc,decreasing = T),]
    return(tt)
}
#rnk
write.table(rnkm(res$condition_old_vs_young),'agingSkinTissue.rnk',quote = F,row.names = F,col.names = F,sep='\t')
