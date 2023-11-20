library(ggplot2)
library(viridis)
library(phateR)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggsci)

Tcell=readRDS("Tcell_in_ann_new_1224.rds")

###fig 3A####
DimPlot(Tcell, reduction = "umap",label = F, group.by = "cluster_name")+scale_color_igv()
ggsave("fig 3A.pdf",width = 6,height = 4,units = "in")

###fig 3B####
DimPlot(Tcell,seed = 1, reduction = "umap",label = F,group.by = 'lesions',label.size = 6) +scale_color_igv()
ggsave("fig 3B.pdf",width = 5.5,height = 4,units = "in")


#fig S7A####
metadata=Tcell@meta.data
Tcellcount=table(metadata$orig.ident,metadata$cluster_name)
Tcellsum=apply(Tcellcount,1,sum)
Tcellfrac=apply(Tcellcount,2,function(x,y){
  return(x/y)
},y=Tcellsum)

EMTome_list=clusterProfiler::read.gmt("EMTome_signatures.gmt")
EMTome_list$gene_new=trimws(EMTome_list$gene)
length(unique(EMTome_list$gene_new)) 


library(GSVA)
gs=list() #gs: gene set
gs[[1]]=unique(EMTome_list$gene_new)
Tdata=GetAssayData(Tcell)
Tdata=as.matrix(Tdata)
celltype=colnames(Tcellfrac)
sampleID=rownames(Tcellfrac)
genedata=do.call(cbind,lapply(sampleID, function(s,celltype,metadata,Tdata){
  a=do.call(cbind,lapply(celltype, function(cc,metadata,Tdata){
    cellID=row.names(metadata[metadata$orig.ident==s&metadata$cluster_name==cc,])
    subdata=Tdata[,match(cellID,colnames(Tdata))]
    return(apply(subdata, 1, mean))
  },metadata,Tdata))
  colnames(a)=paste(s,celltype,sep=":")
  return(a)
},celltype,metadata,Tdata))

genedata[1:6,1:6]
gsva.es <- gsva(genedata, gs, verbose=FALSE)  #es:enriment score
EMTgsva=data.frame(EMTscore=as.numeric(gsva.es),sampleID=do.call(rbind,strsplit(colnames(gsva.es),split = ":"))[,1],celltype=do.call(rbind,strsplit(colnames(gsva.es),split = ":"))[,2])
EMTgsva$lession=metadata$lesions[match(EMTgsva$sampleID,metadata$orig.ident)]
EMTgsva$patient=metadata$patient[match(EMTgsva$sampleID,metadata$orig.ident)]

cluster_name = levels(factor(Tcell@meta.data$cluster_name))
EMTgsva$celltype<-factor(EMTgsva$celltype,level=cluster_name)

EMTgsva = EMTgsva%>%mutate(a=as.numeric(celltype),
                           b=ifelse(lession=='Colon_tumor', -0.2, 0.2)) 
EMTgsva$lession=factor(EMTgsva$lession, levels = c('Colon_tumor','Liver_meta'),ordered = T)
EMTgsva$celltype=factor(EMTgsva$celltype, levels = c('CD4_IL17A_Th17','CD4_FOXP3_Treg','CD4_CCR7_Tn-like1',
                                                     'Cycling_T','CD4_IL7R_Tn-like2','CD8_CXCL13_Tex',
                                                     'CD8_GZMK_Tem','CD8_SLC4A10_MAIT','TYROBP_NK&NKT'),ordered = T)

ggplot(EMTgsva,aes(x=celltype,y=EMTscore,color=lession))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_igv()+
  geom_jitter(aes(a+b, EMTscore, color=lession), position=position_jitter(width = 0.1, height=0),alpha=0.7,size=1)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 45,vjust=0.9,hjust=0.9),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black")) + 
  ggpubr::stat_compare_means(method = "wilcox.test", method.args = list(alternative = "greater"),size=2.5,label = "p.format",angle=0,vjust = 0)+   #仅显示P值，不显示方法
  labs(title = "")+xlab('')+ylab('EMT relatedness')
ggsave("fig S7A.pdf",width = 7,height = 4,units = "in")



