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
fpkm=getfpkm('AP1_inhibitor/RNA/')
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
