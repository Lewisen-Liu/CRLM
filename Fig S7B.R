library(Seurat)
library(dplyr)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(GSVA)

Tcell=readRDS("Tcell_in_ann_new_1224.rds")

Regev.CRC.intratumor.immune=readRDS("Regev_CRC_intratumor_immune.rds")
Regev.CRC.intratumor.Tcell=Regev.CRC.intratumor.immune %>% subset(clTopLevel=='TNKILC') 
Regev.CRC.intratumor.Th17=Regev.CRC.intratumor.Tcell %>% subset(cl295v11SubShort=='cTNI05')  
Regev.CRC.intratumor.Treg=Regev.CRC.intratumor.Tcell %>% subset(cl295v11SubShort=='cTNI08')  

Regev.CRC.adj.immune=readRDS("Regev_CRC_adjnormal_immune.rds")
Regev.CRC.adj.Tcell=Regev.CRC.adj.immune %>% subset(clTopLevel=='TNKILC')
Regev.CRC.adj.Th17=Regev.CRC.adj.Tcell %>% subset(cl295v11SubShort=='cTNI05')  
Regev.CRC.adj.Treg=Regev.CRC.adj.Tcell %>% subset(cl295v11SubShort=='cTNI08') 

zzm_HCC_PT=readRDS("zzm_HCC_intratumor_immune.rds")
zzm_HCC_PT_Tcell=subset(zzm_HCC_PT,celltype_global%in%c('Lymphoid-NK','Lymphoid-T','Lymphoid-T-NK-Cycling','ILCs'))
zzm_HCC_PT_Tcell.Treg=zzm_HCC_PT_Tcell %>% subset(celltype_sub=='CD4-C7-FOXP3')

zzm_HCC_PT_Tcell.Th17=zzm_HCC_PT_Tcell%>%subset(IL17A>0)
table(zzm_HCC_PT_Tcell.Th17$Sample) 

Tcell_out=readRDS("Tcell_out.rds")
Tcell_LM=Tcell_out %>% subset(lesions=="Liver_adjacent") 

length(grep("^IG[HLK]",rownames(Tcell_LM)))
IG_name<-rownames(Tcell_LM)[grep("^IG[HLK]",rownames(Tcell_LM))]
IG_name<-c(IG_name,"JCHAIN")
Tcell_LM<-Tcell_LM[!rownames(Tcell_LM) %in% IG_name]

excludeGene<-read.table("excludeGene.txt",header=TRUE)
Tcell_LM<-Tcell_LM[!rownames(Tcell_LM) %in% excludeGene$Gene]


set.seed(1)
Tcell_LM<-NormalizeData(Tcell_LM,normalization.method = "LogNormalize",scale.factor = 10000)
Tcell_LM<-FindVariableFeatures(Tcell_LM, selection.method = "vst", nfeatures = 2000)
Tcell_LM<-ScaleData(Tcell_LM,features = VariableFeatures(Tcell_LM))
Tcell_LM<-RunPCA(Tcell_LM, features = VariableFeatures(object = Tcell_LM),npcs = 50, verbose = TRUE)
ElbowPlot(Tcell_LM,ndims = 50)
Tcell_LM<-FindNeighbors(Tcell_LM, reduction = "pca",dims = 1:30)
Tcell_LM<-FindClusters(Tcell_LM, resolution = 2) 
Tcell_LM<-RunUMAP(Tcell_LM,dims = 1:30,reduction = "pca")
FeaturePlot(Tcell_LM, features=c('CD3D','CD3E','CD4','FOXP3'),reduction = "umap",cols = c("#DFDFDF","#FF0000"))
FeaturePlot(Tcell_LM, features=c('CD3D','CD3E','CD4','IL17A'),reduction = "umap",cols = c("#DFDFDF","#FF0000"))
DimPlot(Tcell_LM,seed = 1, reduction = "umap",label = T,label.size = 6)+scale_color_igv()


CRLM.adj.Treg.LM=Tcell_LM %>% subset(seurat_clusters=='16') 
CRLM.adj.Th17.LM=Tcell_LM%>%subset(IL17A>0)


Tcell_PT=Tcell_out %>% subset(lesions=="Colon_adjacent")
Tcell_PT=Tcell_PT%>%subset(orig.ident%in%c('O5225-3','P5215-3','Q5215-4','R5218-4'))
length(grep("^IG[HLK]",rownames(Tcell_PT)))
IG_name<-rownames(Tcell_PT)[grep("^IG[HLK]",rownames(Tcell_PT))]
IG_name<-c(IG_name,"JCHAIN")
Tcell_PT<-Tcell_PT[!rownames(Tcell_PT) %in% IG_name]
Tcell_PT<-Tcell_PT[!rownames(Tcell_PT) %in% excludeGene$Gene]

set.seed(1)
Tcell_PT<-NormalizeData(Tcell_PT,normalization.method = "LogNormalize",scale.factor = 10000)
Tcell_PT<-FindVariableFeatures(Tcell_PT, selection.method = "vst", nfeatures = 2000)
Tcell_PT<-ScaleData(Tcell_PT,features = VariableFeatures(Tcell_PT))
Tcell_PT<-RunPCA(Tcell_PT, features = VariableFeatures(object = Tcell_PT),npcs = 50, verbose = TRUE)
ElbowPlot(Tcell_PT,ndims = 50)
Tcell_PT<-FindNeighbors(Tcell_PT, reduction = "pca",dims = 1:30)
Tcell_PT<-FindClusters(Tcell_PT, resolution = 2) 
Tcell_PT<-RunUMAP(Tcell_PT,dims = 1:30,reduction = "pca")
FeaturePlot(Tcell_PT, features=c('CD3D','CD3E','CD4','FOXP3'),reduction = "umap",cols = c("#DFDFDF","#FF0000"))
FeaturePlot(Tcell_PT, features=c('CD3D','CD3E','CD4','IL17A'),reduction = "umap",cols = c("#DFDFDF","#FF0000"))

DimPlot(Tcell_PT,seed = 1, reduction = "umap",label = T,label.size = 6)+scale_color_igv()

CRLM.adj.Treg.PT=Tcell_PT %>% subset(seurat_clusters=='14') 
CRLM.adj.Th17.PT=Tcell_PT%>%subset(IL17A>0)

normal_liver=readRDS("normal_liver.rds")
Liver_FOXP3_Treg=normal_liver %>% subset(FOXP3>0)
Liver_IL17A_Th17=normal_liver%>%subset(IL17A>0)




EMTome_list=clusterProfiler::read.gmt("EMTome_signatures.gmt")
EMTome_list$gene_new=trimws(EMTome_list$gene)
length(unique(EMTome_list$gene_new)) 
gs=list() #gs: gene set
gs[[1]]=unique(EMTome_list$gene_new)


EMTscore=function(Seuratdata,samplelist){
  gsvadata<-as.data.frame(matrix(nrow = dim(Seuratdata)[1],ncol = 1))
  rownames(gsvadata)=rownames(Seuratdata)
  for (i in samplelist) {
    tmp.seurat=Seuratdata %>% subset(PID_new%in%i)  
    genedata=as.matrix(GetAssayData(tmp.seurat))
    genedata_mean=as.matrix(rowMeans(genedata))
    colnames(genedata_mean)=i
    gsvadata=cbind(gsvadata,genedata_mean)
  }
  gsvadata=gsvadata[-1]
  gsvadata=as.matrix(gsvadata)
  gsva.es <- gsva(gsvadata, gs, verbose=FALSE) 
}



Tcell.PT=Tcell %>% subset(lesions=='Colon_tumor')
Tcell.LM=Tcell %>% subset(lesions=='Liver_meta') 


head(Tcell.PT)
CD4_CCR7_Tn_like1.PT=Tcell.PT %>% subset(cluster_name=='CD4_CCR7_Tn-like1')  
CD4_CCR7_Tn_like1.PT$PID_new=paste(CD4_CCR7_Tn_like1.PT$orig.ident, '_','CD4_CCR7_Tn_like1.PT')
length(unique(CD4_CCR7_Tn_like1.PT$PID_new))

CD4_FOXP3_Treg.PT=Tcell.PT %>% subset(cluster_name=='CD4_FOXP3_Treg')  
CD4_FOXP3_Treg.PT$PID_new=paste(CD4_FOXP3_Treg.PT$orig.ident, '_','CD4_FOXP3_Treg.PT')
length(unique(CD4_FOXP3_Treg.PT$PID_new))

CD4_IL17A_Th17.PT=Tcell.PT %>% subset(cluster_name=='CD4_IL17A_Th17')  
CD4_IL17A_Th17.PT$PID_new=paste(CD4_IL17A_Th17.PT$orig.ident, '_','CD4_IL17A_Th17.PT')
length(unique(CD4_IL17A_Th17.PT$PID_new))

CD4_IL7R_Tn_like2.PT=Tcell.PT %>% subset(cluster_name=='CD4_IL7R_Tn-like2') 
CD4_IL7R_Tn_like2.PT$PID_new=paste(CD4_IL7R_Tn_like2.PT$orig.ident, '_','CD4_IL7R_Tn_like2.PT')
length(unique(CD4_IL7R_Tn_like2.PT$PID_new))

CD8_CXCL13_Tex.PT=Tcell.PT %>% subset(cluster_name=='CD8_CXCL13_Tex')  
CD8_CXCL13_Tex.PT$PID_new=paste(CD8_CXCL13_Tex.PT$orig.ident, '_','CD8_CXCL13_Tex.PT')
length(unique(CD8_CXCL13_Tex.PT$PID_new))

CD8_GZMK_Tem.PT=Tcell.PT %>% subset(cluster_name=='CD8_GZMK_Tem')  
CD8_GZMK_Tem.PT$PID_new=paste(CD8_GZMK_Tem.PT$orig.ident, '_','CD8_GZMK_Tem.PT')
length(unique(CD8_GZMK_Tem.PT$PID_new))

CD8_SLC4A10_MAIT.PT=Tcell.PT %>% subset(cluster_name=='CD8_SLC4A10_MAIT')  
CD8_SLC4A10_MAIT.PT$PID_new=paste(CD8_SLC4A10_MAIT.PT$orig.ident, '_','CD8_SLC4A10_MAIT.PT')
length(unique(CD8_SLC4A10_MAIT.PT$PID_new))

Cycling_T.PT=Tcell.PT %>% subset(cluster_name=='Cycling_T')  
Cycling_T.PT$PID_new=paste(Cycling_T.PT$orig.ident, '_','Cycling_T.PT')
length(unique(Cycling_T.PT$PID_new))

TYROBP_NK_NKT.PT=Tcell.PT %>% subset(cluster_name=='TYROBP_NK&NKT')  
TYROBP_NK_NKT.PT$PID_new=paste(TYROBP_NK_NKT.PT$orig.ident, '_','TYROBP_NK_NKT.PT')
length(unique(TYROBP_NK_NKT.PT$PID_new))



head(Tcell.LM)
CD4_CCR7_Tn_like1.LM=Tcell.LM %>% subset(cluster_name=='CD4_CCR7_Tn-like1')  
CD4_CCR7_Tn_like1.LM$PID_new=paste(CD4_CCR7_Tn_like1.LM$orig.ident, '_','CD4_CCR7_Tn_like1.LM')
length(unique(CD4_CCR7_Tn_like1.LM$PID_new))

CD4_FOXP3_Treg.LM=Tcell.LM %>% subset(cluster_name=='CD4_FOXP3_Treg')  
CD4_FOXP3_Treg.LM$PID_new=paste(CD4_FOXP3_Treg.LM$orig.ident, '_','CD4_FOXP3_Treg.LM')
length(unique(CD4_FOXP3_Treg.LM$PID_new))

CD4_IL17A_Th17.LM=Tcell.LM %>% subset(cluster_name=='CD4_IL17A_Th17')  
CD4_IL17A_Th17.LM$PID_new=paste(CD4_IL17A_Th17.LM$orig.ident, '_','CD4_IL17A_Th17.LM')
length(unique(CD4_IL17A_Th17.LM$PID_new))

CD4_IL7R_Tn_like2.LM=Tcell.LM %>% subset(cluster_name=='CD4_IL7R_Tn-like2')  
CD4_IL7R_Tn_like2.LM$PID_new=paste(CD4_IL7R_Tn_like2.LM$orig.ident, '_','CD4_IL7R_Tn_like2.LM')
length(unique(CD4_IL7R_Tn_like2.LM$PID_new))

CD8_CXCL13_Tex.LM=Tcell.LM %>% subset(cluster_name=='CD8_CXCL13_Tex')  
CD8_CXCL13_Tex.LM$PID_new=paste(CD8_CXCL13_Tex.LM$orig.ident, '_','CD8_CXCL13_Tex.LM')
length(unique(CD8_CXCL13_Tex.LM$PID_new))

CD8_GZMK_Tem.LM=Tcell.LM %>% subset(cluster_name=='CD8_GZMK_Tem')  
CD8_GZMK_Tem.LM$PID_new=paste(CD8_GZMK_Tem.LM$orig.ident, '_','CD8_GZMK_Tem.LM')
length(unique(CD8_GZMK_Tem.LM$PID_new))

CD8_SLC4A10_MAIT.LM=Tcell.LM %>% subset(cluster_name=='CD8_SLC4A10_MAIT')  
CD8_SLC4A10_MAIT.LM$PID_new=paste(CD8_SLC4A10_MAIT.LM$orig.ident, '_','CD8_SLC4A10_MAIT.LM')
length(unique(CD8_SLC4A10_MAIT.LM$PID_new))

Cycling_T.LM=Tcell.LM %>% subset(cluster_name=='Cycling_T')  
Cycling_T.LM$PID_new=paste(Cycling_T.LM$orig.ident, '_','Cycling_T.LM')
length(unique(Cycling_T.LM$PID_new))

TYROBP_NK_NKT.LM=Tcell.LM %>% subset(cluster_name=='TYROBP_NK&NKT')  
TYROBP_NK_NKT.LM$PID_new=paste(TYROBP_NK_NKT.LM$orig.ident, '_','TYROBP_NK_NKT.LM')
length(unique(TYROBP_NK_NKT.LM$PID_new))


Regev.CRC.intratumor.Th17$PID_new=paste(Regev.CRC.intratumor.Th17$PID,'_','Regev.CRC.intratumor.Th17')

Regev.CRC.intratumor.Treg$PID_new=paste(Regev.CRC.intratumor.Treg$PID,'_','Regev.CRC.intratumor.Treg')

Regev.CRC.adj.Th17$PID_new=paste(Regev.CRC.adj.Th17$PID,'_','Regev.CRC.adj.Th17')

Regev.CRC.adj.Treg$PID_new=paste(Regev.CRC.adj.Treg$PID,'_','Regev.CRC.adj.Treg')

zzm_HCC_PT_Tcell.Treg$PID_new=paste(zzm_HCC_PT_Tcell.Treg$Sample,'_','zzm_HCC_PT_Tcell.Treg') 

zzm_HCC_PT_Tcell.Th17$PID_new=paste(zzm_HCC_PT_Tcell.Th17$Sample,'_','zzm_HCC_PT_Tcell.Th17') 

Liver_FOXP3_Treg$PID_new=paste(Liver_FOXP3_Treg$orig.ident,'_','Liver_FOXP3_Treg') 

Liver_IL17A_Th17$PID_new=paste(Liver_IL17A_Th17$orig.ident,'_','Liver_IL17A_Th17') 

CRLM.adj.Treg.LM$PID_new=paste(CRLM.adj.Treg.LM$orig.ident,'_','CRLM.adj.Treg.LM') 

CRLM.adj.Treg.PT$PID_new=paste(CRLM.adj.Treg.PT$orig.ident,'_','CRLM.adj.Treg.PT') 

CRLM.adj.Th17.LM$PID_new=paste(CRLM.adj.Th17.LM$orig.ident,'_','CRLM.adj.Th17.LM') 

CRLM.adj.Th17.PT$PID_new=paste(CRLM.adj.Th17.PT$orig.ident,'_','CRLM.adj.Th17.PT') 


Tcell_merge=CD4_CCR7_Tn_like1.PT %>% merge(CD4_FOXP3_Treg.PT) %>% merge(CD4_IL17A_Th17.PT) %>% merge(CD4_IL7R_Tn_like2.PT) %>% merge(CD8_CXCL13_Tex.PT) %>% 
  merge(CD8_GZMK_Tem.PT) %>% merge(CD8_SLC4A10_MAIT.PT) %>% merge(Cycling_T.PT) %>% merge(TYROBP_NK_NKT.PT) %>% 
  merge(CD4_CCR7_Tn_like1.LM) %>% merge(CD4_FOXP3_Treg.LM) %>% merge(CD4_IL17A_Th17.LM) %>% merge(CD4_IL7R_Tn_like2.LM) %>% merge(CD8_CXCL13_Tex.LM) %>% 
  merge(CD8_GZMK_Tem.LM) %>% merge(CD8_SLC4A10_MAIT.LM) %>% merge(Cycling_T.LM) %>% merge(TYROBP_NK_NKT.LM) %>% 
  merge(Regev.CRC.intratumor.Th17) %>% merge(Regev.CRC.intratumor.Treg) %>% merge(Regev.CRC.adj.Th17) %>% 
  merge(Regev.CRC.adj.Treg) %>% merge(zzm_HCC_PT_Tcell.Treg) %>% merge(Liver_FOXP3_Treg)%>%
  merge(zzm_HCC_PT_Tcell.Th17)%>%merge(Liver_IL17A_Th17)%>%
  merge(CRLM.adj.Treg.LM)%>%merge(CRLM.adj.Treg.PT)%>%merge(CRLM.adj.Th17.LM)%>%merge(CRLM.adj.Th17.PT)

samplelist=unique(Tcell_merge$PID_new)

Tcell_merge.score=EMTscore(Tcell_merge, samplelist)

Tcell_EMTscore_data=data.frame(
  source=c(rep('CD4_CCR7_Tn_like1.PT',length(unique(CD4_CCR7_Tn_like1.PT$PID_new))), #16
           rep('CD4_FOXP3_Treg.PT',length(unique(CD4_FOXP3_Treg.PT$PID_new))),   #16
           rep('CD4_IL17A_Th17.PT',length(unique(CD4_IL17A_Th17.PT$PID_new))),   #16
           rep('CD4_IL7R_Tn_like2.PT',length(unique(CD4_IL7R_Tn_like2.PT$PID_new))), #16
           rep('CD8_CXCL13_Tex.PT',length(unique(CD8_CXCL13_Tex.PT$PID_new))),   #16
           rep('CD8_GZMK_Tem.PT',length(unique(CD8_GZMK_Tem.PT$PID_new))),   #16
           rep('CD8_SLC4A10_MAIT.PT',length(unique(CD8_SLC4A10_MAIT.PT$PID_new))), #16
           rep('Cycling_T.PT',length(unique(Cycling_T.PT$PID_new))),   #16
           rep('TYROBP_NK_NKT.PT',length(unique(TYROBP_NK_NKT.PT$PID_new))),   #16
           
           rep('CD4_CCR7_Tn_like1.LM',length(unique(CD4_CCR7_Tn_like1.LM$PID_new))), #17
           rep('CD4_FOXP3_Treg.LM',length(unique(CD4_FOXP3_Treg.LM$PID_new))),   #17
           rep('CD4_IL17A_Th17.LM',length(unique(CD4_IL17A_Th17.LM$PID_new))),   #17
           rep('CD4_IL7R_Tn_like2.LM',length(unique(CD4_IL7R_Tn_like2.LM$PID_new))), #17
           rep('CD8_CXCL13_Tex.LM',length(unique(CD8_CXCL13_Tex.LM$PID_new))),   #17
           rep('CD8_GZMK_Tem.LM',length(unique(CD8_GZMK_Tem.LM$PID_new))),   #17
           rep('CD8_SLC4A10_MAIT.LM',length(unique(CD8_SLC4A10_MAIT.LM$PID_new))), #17
           rep('Cycling_T.LM',length(unique(Cycling_T.LM$PID_new))),   #17
           rep('TYROBP_NK_NKT.LM',length(unique(TYROBP_NK_NKT.LM$PID_new))),   #17
           
           rep('Regev.CRC.intratumor.Th17',length(unique(Regev.CRC.intratumor.Th17$PID_new))),  #26
           rep('Regev.CRC.intratumor.Treg',length(unique(Regev.CRC.intratumor.Treg$PID_new))),  #27
           rep('Regev.CRC.adj.Th17',length(unique(Regev.CRC.adj.Th17$PID_new))),   #12
           rep('Regev.CRC.adj.Treg',length(unique(Regev.CRC.adj.Treg$PID_new))),  #14
           rep('zzm_HCC_PT_Tcell.Treg',length(unique(zzm_HCC_PT_Tcell.Treg$PID_new))),  #8
           rep('Liver_FOXP3_Treg',length(unique(Liver_FOXP3_Treg$PID_new))), #9
           rep('zzm_HCC_PT_Tcell.Th17',length(unique(zzm_HCC_PT_Tcell.Th17$PID_new))), #4
           rep('Liver_IL17A_Th17',length(unique(Liver_IL17A_Th17$PID_new))), #7
           
           rep('CRLM.adj.Treg.LM',length(unique(CRLM.adj.Treg.LM$PID_new))), #2
           rep('CRLM.adj.Treg.PT',length(unique(CRLM.adj.Treg.PT$PID_new))), #4
           rep('CRLM.adj.Th17.LM',length(unique(CRLM.adj.Th17.LM$PID_new))), #1
           rep('CRLM.adj.Th17.PT',length(unique(CRLM.adj.Th17.PT$PID_new)))), #3
  
  EMTscore=c(unname(Tcell_merge.score[,1:16]),
             unname(Tcell_merge.score[,17:32]),
             unname(Tcell_merge.score[,33:48]),
             unname(Tcell_merge.score[,49:64]),
             unname(Tcell_merge.score[,65:80]),
             unname(Tcell_merge.score[,81:96]),
             unname(Tcell_merge.score[,97:112]),
             unname(Tcell_merge.score[,113:128]),
             unname(Tcell_merge.score[,129:144]),
             
             unname(Tcell_merge.score[,145:161]),
             unname(Tcell_merge.score[,162:178]),
             unname(Tcell_merge.score[,179:195]),
             unname(Tcell_merge.score[,196:212]),
             unname(Tcell_merge.score[,213:229]),
             unname(Tcell_merge.score[,230:246]),
             unname(Tcell_merge.score[,247:263]),
             unname(Tcell_merge.score[,264:280]),
             unname(Tcell_merge.score[,281:297]), 
             
             unname(Tcell_merge.score[,298:323]),
             unname(Tcell_merge.score[,324:350]),
             unname(Tcell_merge.score[,351:362]),
             unname(Tcell_merge.score[,363:376]),
             unname(Tcell_merge.score[,377:384]),
             unname(Tcell_merge.score[,385:393]),
             
             unname(Tcell_merge.score[,394:397]),
             unname(Tcell_merge.score[,398:404]),
             
             unname(Tcell_merge.score[,405:406]),
             unname(Tcell_merge.score[,407:410]),
             unname(Tcell_merge.score[,411]),
             unname(Tcell_merge.score[,412:414])
  )
)


compaired=list(c('CD4_FOXP3_Treg.PT','CD4_FOXP3_Treg.LM'),c('CD4_IL17A_Th17.PT','CD4_IL17A_Th17.LM'))





EMTgsva=readRDS("/Users/a1502/Desktop/TMP/Tcell_EMTgsva.rds") #输出数据再合并展示
head(EMTgsva)
EMTgsva %>% filter(celltype=='CD4_FOXP3_Treg',lession=='Colon_tumor')

colnames(EMTgsva)[3]='source'

#####Fig S7B_1#######
unique(Tcell_EMTscore_data$source)
tmp1=Tcell_EMTscore_data %>% filter(source%in%c('Regev.CRC.intratumor.Th17','Regev.CRC.adj.Th17','CRLM.adj.Th17.PT'))

tmp2=EMTgsva %>% filter(lession=='Colon_tumor',source=='CD4_IL17A_Th17')
tmp2[,8]=tmp2[,1]
tmp2=tmp2[,c(3,8)]
colnames(tmp2)[2]='EMTscore'

tmp3=rbind(tmp1,tmp2)

tmp3$source[tmp3$source=='CRLM.adj.Th17.PT']='Regev.CRC.adj.Th17'
compaired=list(c('CD4_IL17A_Th17','Regev.CRC.adj.Th17'),c('CD4_IL17A_Th17','Regev.CRC.intratumor.Th17'))

tmp3$source=factor(tmp3$source, levels = c('Regev.CRC.adj.Th17','Regev.CRC.intratumor.Th17','CD4_IL17A_Th17'),ordered = T)
ggplot(tmp3,aes(x=source,y=EMTscore,color=source))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  #scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(colour = "black"),
        axis.text.x = element_text(angle = 30,vjust=0.5,hjust=0.5), #angle = 45,vjust=0.9,hjust=0.9
        plot.title = element_text(hjust = 0.5,size = 12,colour = "black"))+
  labs(title = "Th17 EMT relatedness (Colon)")+ylab('EMT relatedness')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('ANT','nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('ANT','nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))+xlab('')
ggsave("Fig S7B_1.pdf",width = 4,height = 6,units = "in")


#####Fig S7B_3#######
unique(Tcell_EMTscore_data$source)
tmp1=Tcell_EMTscore_data %>% filter(source%in%c('zzm_HCC_PT_Tcell.Th17','Liver_IL17A_Th17','CRLM.adj.Th17.LM'))

tmp2=EMTgsva %>% filter(lession=='Liver_meta',source=='CD4_IL17A_Th17')
tmp2[,8]=tmp2[,1]
tmp2=tmp2[,c(3,8)]
colnames(tmp2)[2]='EMTscore'

tmp3=rbind(tmp1,tmp2)

tmp3$source[tmp3$source=='CRLM.adj.Th17.LM']='Liver_IL17A_Th17'
compaired=list(c('CD4_IL17A_Th17','Liver_IL17A_Th17'),c('CD4_IL17A_Th17','zzm_HCC_PT_Tcell.Th17'))

tmp3$source=factor(tmp3$source, levels = c('Liver_IL17A_Th17','zzm_HCC_PT_Tcell.Th17','CD4_IL17A_Th17'),ordered = T)
ggplot(tmp3,aes(x=source,y=EMTscore,color=source))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  #scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(colour = "black"),
        axis.text.x = element_text(angle = 30,vjust=0.5,hjust=0.5), #angle = 45,vjust=0.9,hjust=0.9
        plot.title = element_text(hjust = 0.5,size = 12,colour = "black"))+
  labs(title = "Th17 EMT relatedness (Liver)")+ylab('EMT relatedness')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('ANT','HCC','CRLM'))+
  scale_colour_discrete(labels=c('ANT','HCC','CRLM'),type = c('#1f77b4','#ff7f0e','#2ca02c'))+xlab('')

ggsave("Fig S7B_3.pdf",width = 4,height = 6,units = "in")



#####Fig S7B_2#######
unique(Tcell_EMTscore_data$source)
tmp1=Tcell_EMTscore_data %>% filter(source%in%c('Regev.CRC.intratumor.Treg','Regev.CRC.adj.Treg','CRLM.adj.Treg.PT'))

tmp2=EMTgsva %>% filter(lession=='Colon_tumor',source=='CD4_FOXP3_Treg')
tmp2[,8]=tmp2[,1]
tmp2=tmp2[,c(3,8)]
colnames(tmp2)[2]='EMTscore'

tmp3=rbind(tmp1,tmp2)

tmp3$source[tmp3$source=='CRLM.adj.Treg.PT']='Regev.CRC.adj.Treg'
compaired=list(c('CD4_FOXP3_Treg','Regev.CRC.adj.Treg'),c('CD4_FOXP3_Treg','Regev.CRC.intratumor.Treg'))

tmp3$source=factor(tmp3$source, levels = c('Regev.CRC.adj.Treg','Regev.CRC.intratumor.Treg','CD4_FOXP3_Treg'),ordered = T)
ggplot(tmp3,aes(x=source,y=EMTscore,color=source))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(colour = "black"),
        axis.text.x = element_text(angle = 30,vjust=0.5,hjust=0.5), #angle = 45, vjust=0.9, hjust=0.9
        plot.title = element_text(hjust = 0.5,size = 12, colour = "black"))+
  labs(title = "Treg EMT relatedness (Colon)")+ylab('EMT relatedness')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('ANT','nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('ANT','nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))+xlab('')
ggsave("Fig S7B_2.pdf",width = 4,height = 6,units = "in")




#####Fig S7B_4#######
unique(Tcell_EMTscore_data$source)
tmp1=Tcell_EMTscore_data %>% filter(source%in%c('zzm_HCC_PT_Tcell.Treg','Liver_FOXP3_Treg','CRLM.adj.Treg.LM'))

tmp2=EMTgsva %>% filter(lession=='Liver_meta',source=='CD4_FOXP3_Treg')
tmp2[,8]=tmp2[,1]
tmp2=tmp2[,c(3,8)]
colnames(tmp2)[2]='EMTscore'

tmp3=rbind(tmp1,tmp2)

tmp3$source[tmp3$source=='CRLM.adj.Treg.LM']='Liver_FOXP3_Treg'
compaired=list(c('CD4_FOXP3_Treg','Liver_FOXP3_Treg'),c('CD4_FOXP3_Treg','zzm_HCC_PT_Tcell.Treg'))

tmp3$source=factor(tmp3$source, levels = c('Liver_FOXP3_Treg','zzm_HCC_PT_Tcell.Treg','CD4_FOXP3_Treg'),ordered = T)
ggplot(tmp3,aes(x=source,y=EMTscore,color=source))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(colour = "black"),
        axis.text.x = element_text(angle = 30,vjust=0.5,hjust=0.5), #angle = 45,vjust=0.9,hjust=0.9
        plot.title = element_text(hjust = 0.5,size = 12,colour = "black"))+
  labs(title = "Treg EMT relatedness (Liver)")+ylab('EMT relatedness')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('ANT','HCC','mCRC'))+
  scale_colour_discrete(labels=c('ANT','HCC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))+xlab('')
ggsave("Fig S7B_4.pdf",width = 4,height = 6,units = "in")

