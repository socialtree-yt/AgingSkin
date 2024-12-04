depp<-c("Seurat",'ggplot2')
sapply(depp, library, character.only = TRUE)
scRNA.integrated=readRDS('../05.1-subcluster/subtype.harmony.rds')
options(repr.plot.width=8, repr.plot.height=6)
colors<-c("#4d648d", "#0450fb", "#11aac4", "#42e8f3", "#AEC7E8", "#2CA02C", "#98DF8A", "#9eccaf", "#daf400", "#983b59", "#e81f3f", "#ff8b94", "#ffd3b6", "#f9ae34", "#ffdb00", "#723584", "#9264eb", "#ff00ff", "#E377C2", "#de94e4", "#F7B6D2", "#C5B0D5", "#8C564B", "#C49C94", "#BCBD22", "#DBDB8D", "#7F7F7F", "#C7C7C7", "#a7a7b7", "#9999FF",  "#CC6633", "#990066", "#003333", "#996666","#CF6633")
DimPlot(scRNA.integrated, reduction = "umap", label = TRUE, repel = FALSE, cols=colors[1:length(unique(scRNA.integrated@active.ident))]) + NoLegend() + ggtitle("UMAP plot: colored by cluster")
DimPlot(scRNA.integrated, reduction = "umap",group.by = 'group', label = TRUE, repel = FALSE, cols=colors[1:length(unique(scRNA.integrated@active.ident))]) + NoLegend() + ggtitle("UMAP plot: colored by cluster")
DimPlot(scRNA.integrated, reduction = "umap",group.by = 'sample', label = TRUE, repel = FALSE, cols=colors[1:length(unique(scRNA.integrated@active.ident))]) + NoLegend() + ggtitle("UMAP plot: colored by cluster")
tt=list(
    c('DSC1','KRT2','IVL','TGM3','FLG'), # GC   granular cell
    c('KRT10','KRT1','DSG1','CDH1'),      # SC  spinous cell   also expressed in GC
    c('KRT14','KRT5','CDH3','TP63'),   #BC     basal cell   
    c('KRT15','ITGB1') ##  Epsc  epidermal stem cell  'KRT5',
    )    
names(tt)=c('GC',"SC",'BC','KC')
options(repr.plot.width=16, repr.plot.height=10)
DotPlot(scRNA.integrated, features = tt, group.by = "seurat_clusters")+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("") #+
#theme(legend.position = 'right',axis.title = element_blank(),legend.text=element_blank(),axis.text = element_blank(),
#      axis.text.x = element_blank(),legend.title=element_blank(),legend.key.height=unit(2,"line"))+scale_size(range=c(0, 9))
ggsave(filename = "marker.png",device = "png",width = 64,height = 33,units = "cm")
#DimPlot(object = scRNA.integrated, reduction = "umap", label = TRUE, repel = TRUE, split.by = "group2")
# 16,15,14,13 ,12,11 ,10 ,9? ,8 ,7? ,6,5,4?, 3, 2, 1, 0.
## Assigning cell type identity to clusters
#6,12
#```{r cellType_annotation, message=FALSE,warning=FALSE, echo=FALSE, results='hide', fig.height=6, fig.width=12}
#new.cluster.ids<-c("0: Cycling Epithelial Cells", "1: Surface Ectoderm Cells", "2: Posterior Placodal Ectoderm",  "3: Cycling Epithelial Cells",  "4: Surface Ectoderm Cells", "5: Anterior Placode Ectoderm",  "6: Posterior Placode Ectoderm", "7: Surface Ectoderm Cells", "8: Neuroectoderm", "9: CNCC Mesenchyme", "10: Surface Ectoderm Cells", "11: Cycling Epithelial Cells", "12: Surface Ectoderm Cells", "13: CNCC Mesenchyme", "14: Surface Ectoderm", "15: Cycling Epithelial Cells", "16: Surface Ectoderm Cells", "17: CNCCs", "18: Posterior Placodal Ectoderm", "19: Cycling Neuroectoderm", "20: Low Mito Cells", "21: Neuroectoderm", "22: Neuroblasts", "23: Neuroblasts")
new.cluster.ids<-c('0:BC','1:BC','2:SC','3:SC','4:SC','5:SC','6:SC_BC',
                   '7:BC','8:SC_BC','9:SC','10:SC','11:BC','12:GC_SC','13:SC_BC',
                   '14:BC','15:SC','16:SC_BC','17:SC','18:SC','19:BC','20:BC','21:BC','22:SC','23:SC_BC','24:GC_SC',
                   '25:SC','26:BC','27:BC','28:SC','29:SC_BC','30:SC','31:SC_BC','32:BC','33:SC_BC',
                   '34:BC','35:GC_SC','36:SC_BC'
                  )

scRNA.temp<-scRNA.integrated
names(new.cluster.ids) <- levels(scRNA.temp)
scRNA.temp <- RenameIdents(scRNA.temp, new.cluster.ids)
scRNA.temp[['celltype']]=Idents(scRNA.temp)
#DimPlot(scRNA.temp, reduction = "umap", label = TRUE, cols=colors[1:length(unique(scRNA.temp@active.ident))])
#ggsave(filename = "cell_type.png",device = "png",width = 34,height = 28,units = "cm")
scRNA.integrated=AddModuleScore(scRNA.integrated,features = list(tt[[1]]),name=names(tt)[1])
scRNA.integrated=AddModuleScore(scRNA.integrated,features = list(tt[[2]]),name=names(tt)[2])
scRNA.integrated=AddModuleScore(scRNA.integrated,features = list(tt[[3]]),name=names(tt)[3])
options(repr.plot.width=12, repr.plot.height=12)
VlnPlot(scRNA.integrated,features = 'GC1')
VlnPlot(scRNA.integrated,features = 'SC1')
VlnPlot(scRNA.integrated,features = 'BC1')
options(repr.plot.width=8, repr.plot.height=8)
p1=FeaturePlot(scRNA.integrated,features = 'GC1',min.cutoff = 0,cols = c('#FFFFCC','blue'))+theme(legend.position = 'right',legend.text=element_blank(),legend.key.height=unit(3,"line"),axis.title = element_blank(),
      axis.text = element_blank()
    )
ggsave(filename = "feature.celltype.1.png",p1,device = "png",width = 13.5,height = 15,units = "cm")
p2=FeaturePlot(scRNA.integrated,features = 'SC1',min.cutoff = 0,cols = c('#FFFFCC','blue'))+theme(legend.position = 'right',legend.text=element_blank(),legend.key.height=unit(3,"line"),axis.title = element_blank(),
      axis.text = element_blank()
    )
ggsave(filename = "feature.celltype.2.png",p2,device = "png",width = 13.5,height = 15,units = "cm")
p3=FeaturePlot(scRNA.integrated,features = 'BC1',min.cutoff = 0,cols = c('#FFFFCC','blue'))+theme(legend.position = 'right',legend.text=element_blank(),legend.key.height=unit(3,"line"),axis.title = element_blank(),
      axis.text = element_blank()
    )
ggsave(filename = "feature.celltype.3.png",p3,device = "png",width = 13.5,height = 15,units = "cm")
library(stringr)
ident_id = scRNA.temp[["celltype"]]
newid=str_split_fixed(ident_id[,1],':',n=2)[,2]
scRNA.temp <- AddMetaData(object = scRNA.temp, metadata = newid, col.name = "celltype2") 
Idents(scRNA.temp)=scRNA.temp$celltype2
options(repr.plot.width=8, repr.plot.height=8)
DimPlot(scRNA.temp, reduction = "umap", label = F, cols=colors[1:length(unique(scRNA.temp@active.ident))])+xlim(-10,10)+
        theme(panel.background = element_blank(),axis.line.x = element_line(color="black",size=0.8),axis.line.y = element_line(color="black",size=0.8),
        axis.text.x = element_text(size=10),axis.text.y = element_text(size=10),axis.line = element_line(colour = "black"),
             axis.ticks.x=element_line(color="black",size=0.8),
      axis.ticks.y=element_line(color="black",size=0.8),axis.ticks.length.x = unit(0.2,'cm'), 
      axis.ticks.length.y = unit(0.2,'cm'))
#DimPlot(scRNA.temp, reduction = "umap", label = F,group.by ='seurat_clusters' , cols=colors[1:length(unique(scRNA.integrated@active.ident))])
ggsave(filename = "Figure1B_cell_type.png",device = "png",width = 12,height = 10,units = "cm")
saveRDS(scRNA.temp, file = "scRNA_annodata.KC.rds")
