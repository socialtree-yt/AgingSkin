### merge samples
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
fpkm1=getfpkm('EPSC_GSE92967')
fpkm2=getfpkm('TPA_EPSC_GSE171596/RNA-seq')
fpkm=merge(fpkm1,fpkm2,by = 'Geneid')
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


### calculate degs
fpkm=readRDS('counts.rds')

col=data.frame(row.names = colnames(fpkm),'external_id'=colnames(fpkm),age2=c(rep('D6Ctrl',2),rep('D6IMQl',2),rep('D30Ctrl',2),
                                                                              rep('D30IMQ',2),rep('D30_Ctrl_6hr_TPA',2),rep('D30_IMQ_6hr_TPA',2)))
#head(col,2)
countData=fpkm
colData=col
condition=factor(colData$age2)
#condition = relevel(condition, c("D6Ctrl"))
colData=data.frame(row.names=rownames(colData),condition)

res_total=list()
library(edgeR)
countData <- countData[rowMeans(countData)>1,]
countData=round(countData)
de=function(countData,colData,condition2){
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition2)
    #dds$condition <- relevel(dds$condition, ref = "young")
    head(dds)
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
    return(res)
}

for(lev in levels(condition)){
    condition2 = relevel(condition, c(lev))
    colData2=data.frame(row.names=rownames(colData),condition2)
###################################去除低表达和表达比例低的gene;并取整
    res_total[[lev]]=de(countData,colData2,condition2)
}
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
#dds$condition <- relevel(dds$condition, ref = "young")
########### 计算差异倍数
dds1 <- DESeq(dds) 
normfpkm=counts(dds1, normalized=TRUE)
write.table(normfpkm,'norm.exp.count.txt',row.names = T,quote = F,col.names = T,sep='\t')

### get rnk
rnkm=function(deg){
    tt=data.frame('gene'=rownames(deg),'logfc'=deg$log2FoldChange)
    tt=tt[order(tt$logfc,decreasing = T),]
    return(tt)
}
#rnk
write.table(rnkm(res_total$D30Ctrl$condition2_D30_Ctrl_6hr_TPA_vs_D30Ctrl),'D30_Ctrl_6hr_TPA.D30Ctrl.rnk',quote = F,row.names = F,col.names = F,sep='\t')
write.table(rnkm(res_total$D30_Ctrl_6hr_TP$condition2_D30_IMQ_6hr_TPA_vs_D30_Ctrl_6hr_TPA),'D30_IMQ_6hr_TPA.D30_Ctrl_6hr_TPA.rnk',quote = F,row.names = F,col.names = F,sep='\t')
write.table(rnkm(res_total$D30IMQ$condition2_D30_IMQ_6hr_TPA_vs_D30IMQ),'D30_IMQ_6hr_TPA.D30IMQ.rnk',quote = F,row.names = F,col.names = F,sep='\t')
write.table(rnkm(res_total$D30Ctrl$condition2_D30IMQ_vs_D30Ctrl),'D30IMQ.D30Ctrl.rnk',quote = F,row.names = F,col.names = F,sep='\t')
write.table(rnkm(res_total$D6Ctrl$condition2_D6IMQl_vs_D6Ctrl),'D6IMQ.D6Ctrl.rnk',quote = F,row.names = F,col.names = F,sep='\t')

#degs
cut=1
write.table(rownames(res_total$D30Ctrl$condition2_D30_Ctrl_6hr_TPA_vs_D30Ctrl)[res_total$D30Ctrl$condition2_D30_Ctrl_6hr_TPA_vs_D30Ctrl$pvalue<0.05 & 
                                                                   res_total$D30Ctrl$condition2_D30_Ctrl_6hr_TPA_vs_D30Ctrl$log2FoldChange>=cut]
,'D30_Ctrl_6hr_TPA.D30Ctrl.up.txt',quote = F,row.names = F,col.names = F,sep='\t')
write.table(rownames(res_total$D30Ctrl$condition2_D30_Ctrl_6hr_TPA_vs_D30Ctrl)[res_total$D30Ctrl$condition2_D30_Ctrl_6hr_TPA_vs_D30Ctrl$pvalue<0.05 & 
                                                                   res_total$D30Ctrl$condition2_D30_Ctrl_6hr_TPA_vs_D30Ctrl$log2FoldChange<(-cut)]
,'D30_Ctrl_6hr_TPA.D30Ctrl.down.txt',quote = F,row.names = F,col.names = F,sep='\t')

write.table(rownames(res_total$D30_Ctrl_6hr_TP$condition2_D30_IMQ_6hr_TPA_vs_D30_Ctrl_6hr_TPA)[res_total$D30_Ctrl_6hr_TP$condition2_D30_IMQ_6hr_TPA_vs_D30_Ctrl_6hr_TPA$pvalue<0.05 & 
                                                            res_total$D30_Ctrl_6hr_TP$condition2_D30_IMQ_6hr_TPA_vs_D30_Ctrl_6hr_TPA$log2FoldChange>=cut]
,'D30_IMQ_6hr_TPA.D30_Ctrl_6hr_TPA.up.txt',quote = F,row.names = F,col.names = F,sep='\t')
write.table(rownames(res_total$D30_Ctrl_6hr_TP$condition2_D30_IMQ_6hr_TPA_vs_D30_Ctrl_6hr_TPA)[res_total$D30_Ctrl_6hr_TP$condition2_D30_IMQ_6hr_TPA_vs_D30_Ctrl_6hr_TPA$pvalue<0.05 & 
                                                            res_total$D30_Ctrl_6hr_TP$condition2_D30_IMQ_6hr_TPA_vs_D30_Ctrl_6hr_TPA$log2FoldChange<(-cut)]
,'D30_IMQ_6hr_TPA.D30_Ctrl_6hr_TPA.down.txt',quote = F,row.names = F,col.names = F,sep='\t')

write.table(rownames(res_total$D30IMQ$condition2_D30_IMQ_6hr_TPA_vs_D30IMQ)[res_total$D30IMQ$condition2_D30_IMQ_6hr_TPA_vs_D30IMQ$pvalue<0.05 & 
                                                                   res_total$D30IMQ$condition2_D30_IMQ_6hr_TPA_vs_D30IMQ$log2FoldChange>=cut]
,'D30_IMQ_6hr_TPA.D30IMQ.up.txt',quote = F,row.names = F,col.names = F,sep='\t')
write.table(rownames(res_total$D30IMQ$condition2_D30_IMQ_6hr_TPA_vs_D30IMQ)[res_total$D30IMQ$condition2_D30_IMQ_6hr_TPA_vs_D30IMQ$pvalue<0.05 & 
                                                                   res_total$D30IMQ$condition2_D30_IMQ_6hr_TPA_vs_D30IMQ$log2FoldChange<(-cut)]
,'D30_IMQ_6hr_TPA.D30IMQ.down.txt',quote = F,row.names = F,col.names = F,sep='\t')

write.table(rownames(res_total$D30Ctrl$condition2_D30IMQ_vs_D30Ctrl)[res_total$D30Ctrl$condition2_D30IMQ_vs_D30Ctrl$pvalue<0.05 & 
                                                                   res_total$D30Ctrl$condition2_D30IMQ_vs_D30Ctrl$log2FoldChange>=cut]
,'D30IMQ.D30Ctrl.up.txt',quote = F,row.names = F,col.names = F,sep='\t')
write.table(rownames(res_total$D30Ctrl$condition2_D30IMQ_vs_D30Ctrl)[res_total$D30Ctrl$condition2_D30IMQ_vs_D30Ctrl$pvalue<0.05 & 
                                                                   res_total$D30Ctrl$condition2_D30IMQ_vs_D30Ctrl$log2FoldChange<(-cut)]
,'D30IMQ.D30Ctrl.down.txt',quote = F,row.names = F,col.names = F,sep='\t')

write.table(rownames(res_total$D6Ctrl$condition2_D6IMQl_vs_D6Ctrl)[res_total$D6Ctrl$condition2_D6IMQl_vs_D6Ctrl$pvalue<0.05 & 
                                                                   res_total$D6Ctrl$condition2_D6IMQl_vs_D6Ctrl$log2FoldChange>=cut]
,'D6IMQ.D6Ctrl.up.txt',quote = F,row.names = F,col.names = F,sep='\t')
write.table(rownames(res_total$D6Ctrl$condition2_D6IMQl_vs_D6Ctrl)[res_total$D6Ctrl$condition2_D6IMQl_vs_D6Ctrl$pvalue<0.05 & 
                                                                   res_total$D6Ctrl$condition2_D6IMQl_vs_D6Ctrl$log2FoldChange<(-cut)]
,'D6IMQ.D6Ctrl.down.txt',quote = F,row.names = F,col.names = F,sep='\t')

