library(Seurat)
library(dplyr)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(GSVA)
Regev.CRC.intratumor.immune=readRDS("Regev_CRC_intratumor_immune.rds")
Regev.CRC.adj.immune=readRDS("Regev_CRC_adjnormal_immune.rds")
Regev.CRC.intratumor.epi=readRDS("Regev_CRC_intratumor_tumor.rds")

Tcell_out=readRDS("Tcell_out.rds")
Tcell_LM=Tcell_out %>% subset(lesions=="Liver_adjacent") 
unique(Tcell_LM$orig.ident)

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
table(CRLM.adj.Treg.LM$orig.ident)
CRLM.adj.Th17.LM=Tcell_LM%>%subset(IL17A>0)


Tcell_PT=Tcell_out %>% subset(lesions=="Colon_adjacent")
unique(Tcell_PT$orig.ident)
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
table(CRLM.adj.Treg.PT$orig.ident)
CRLM.adj.Th17.PT=Tcell_PT%>%subset(IL17A>0)



########Fig S7C_1########
Regev.CRC.intratumor.Tcell=Regev.CRC.intratumor.immune %>% subset(clTopLevel=='TNKILC') 
Regev.CRC.adj.Tcell=Regev.CRC.adj.immune %>% subset(clTopLevel=='TNKILC') 

Regev.CRC.intratumor.Th17.prop=table(Regev.CRC.intratumor.Tcell$cl295v11SubShort,Regev.CRC.intratumor.Tcell$PID)['cTNI05',]/table(Regev.CRC.intratumor.Tcell$PID)
Regev.CRC.adj.Th17.prop=table(Regev.CRC.adj.Tcell$cl295v11SubShort,Regev.CRC.adj.Tcell$PID)['cTNI05',]/table(Regev.CRC.adj.Tcell$PID)

CRLM.adj.Th17.PT$orig.ident=factor(CRLM.adj.Th17.PT$orig.ident,levels = unique(factor(Tcell_PT$orig.ident)),ordered = T)
CRLM.adj.Th17.PT.prop=table(CRLM.adj.Th17.PT$orig.ident)/table(Tcell_PT$orig.ident)

Tcell=readRDS("Tcell_in_ann_new_1224.rds")
Tcell.PT=Tcell %>% subset(lesions=='Colon_tumor')
Tcell.PT.Th17.prop=table(Tcell.PT$cluster_name,Tcell.PT$orig.ident)['CD4_IL17A_Th17',]/table(Tcell.PT$orig.ident)
unname(Tcell.PT.Th17.prop)

Th17_data=data.frame(
  source=c(rep('met_CRC',16),rep('nonmet_CRC',27),rep('adj_normal',14+4)),
  proportion=c(unname(Tcell.PT.Th17.prop),unname(Regev.CRC.intratumor.Th17.prop),unname(Regev.CRC.adj.Th17.prop),unname(CRLM.adj.Th17.PT.prop))
)

Th17_data$source<-factor(Th17_data$source,levels = c('adj_normal','nonmet_CRC','met_CRC'), ordered=T) 
ggplot(Th17_data,aes(x=source,y=proportion,color=source))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(colour = "black"),
        axis.text.x = element_text(angle = 30,vjust=0.5,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,colour = "black"))+ #face='bold
  labs(title = "Th17 proportion comparison (PT)")+ylab('Th17 proportion')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('ANT','nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('ANT','nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("Fig S7C_1.pdf",width = 4,height = 6,units = "in")





########Fig S7C_2########
Regev.CRC.intratumor.Treg.prop=table(Regev.CRC.intratumor.Tcell$cl295v11SubShort,Regev.CRC.intratumor.Tcell$PID)['cTNI08',]/table(Regev.CRC.intratumor.Tcell$PID)
Regev.CRC.adj.Treg.prop=table(Regev.CRC.adj.Tcell$cl295v11SubShort,Regev.CRC.adj.Tcell$PID)['cTNI08',]/table(Regev.CRC.adj.Tcell$PID)
Tcell.PT.Treg.prop=table(Tcell.PT$cluster_name,Tcell.PT$orig.ident)['CD4_FOXP3_Treg',]/table(Tcell.PT$orig.ident)

CRLM.adj.Treg.PT$orig.ident=factor(CRLM.adj.Treg.PT$orig.ident,levels = unique(factor(Tcell_PT$orig.ident)),ordered = T)
CRLM.adj.Treg.PT.prop=table(CRLM.adj.Treg.PT$orig.ident)/table(Tcell_PT$orig.ident)

Treg_PT_data=data.frame(
  source=c(rep('met_CRC',16),rep('nonmet_CRC',27),rep('adj_normal',14+4)),
  proportion=c(unname(Tcell.PT.Treg.prop),unname(Regev.CRC.intratumor.Treg.prop),unname(Regev.CRC.adj.Treg.prop),unname(CRLM.adj.Treg.PT.prop))
)

Treg_PT_data$source<-factor(Treg_PT_data$source,levels = c('adj_normal','nonmet_CRC','met_CRC'), ordered=T) 
compaired=list(c('met_CRC','adj_normal'),c('met_CRC','nonmet_CRC'))

ggplot(Treg_PT_data,aes(x=source,y=proportion,color=source))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(colour = "black"),
        axis.text.x = element_text(angle = 30,vjust=0.5,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,colour = "black"))+
  labs(title = "Treg proportion comparison (PT)")+ylab('Treg proportion')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('ANT','nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('ANT','nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("Fig S7C_2.pdf",width = 4,height = 6,units = "in")




########Fig S7C_4########
zzm_HCC_PT=readRDS("zzm_HCC_intratumor_immune.rds") 
zzm_HCC_PT_Tcell=subset(zzm_HCC_PT,celltype_global%in%c('Lymphoid-NK','Lymphoid-T','Lymphoid-T-NK-Cycling','ILCs'))
zzm_HCC_PT_Tcell.Treg.prop=table(zzm_HCC_PT_Tcell$celltype_sub,zzm_HCC_PT_Tcell$Sample)['CD4-C7-FOXP3',]/table(zzm_HCC_PT_Tcell$Sample) 
zzm_HCC_Th17=zzm_HCC_PT_Tcell %>% subset(IL17A> 0)

Tcell_LM=Tcell_out %>% subset(lesions=="Liver_adjacent") 
Tcell_LM<-NormalizeData(Tcell_LM,normalization.method = "LogNormalize",scale.factor = 10000)
Tcell_LM<-FindVariableFeatures(Tcell_LM, selection.method = "vst", nfeatures = 2000)
Tcell_LM<-ScaleData(Tcell_LM,features = VariableFeatures(Tcell_LM))
Tcell_LM<-RunPCA(Tcell_LM, features = VariableFeatures(object = Tcell_LM),npcs = 50, verbose = TRUE)
ElbowPlot(Tcell_LM,ndims = 50)
Tcell_LM<-FindNeighbors(Tcell_LM, reduction = "pca",dims = 1:30)
Tcell_LM<-FindClusters(Tcell_LM, resolution = 1) #C16为Treg
Tcell_LM<-RunUMAP(Tcell_LM,dims = 1:30,reduction = "pca")

head(Tcell_LM)
DimPlot(Tcell_LM,seed = 1, reduction = "umap",label = T,label.size = 6)+scale_color_igv()
FeaturePlot(Tcell_LM, features=c('CD3D','CD3E','CD4','FOXP3'),reduction = "umap",cols = c("#DFDFDF","#FF0000"))
FeaturePlot(Tcell_LM, features=c('FOXP3','IL17A'),reduction = "umap",cols = c("#DFDFDF","#FF0000"))

unique(Tcell_LM$orig.ident)
CRLM_adj.Treg.prop=table(Tcell_LM$seurat_clusters,Tcell_LM$orig.ident)['16',]/table(Tcell_LM$orig.ident) 

Tcell.LM=Tcell %>% subset(lesions=='Liver_meta')
unique(Tcell.LM$orig.ident)
CRLM_tumor.Treg.prop=table(Tcell.LM$cluster_name,Tcell.LM$orig.ident)['CD4_FOXP3_Treg',]/table(Tcell.LM$orig.ident) 

normal_liver=readRDS("normal_liver.rds")
head(normal_liver)
unique(normal_liver$orig.ident) 
unique(normal_liver$batch)

Liver_FOXP3_Treg=normal_liver %>% subset(FOXP3>0)
table(Liver_FOXP3_Treg$orig.ident) 
Liver_FOXP3_Treg$orig.ident=factor(Liver_FOXP3_Treg$orig.ident,levels = unique(factor(normal_liver$orig.ident)),ordered = T)
Liver_FOXP3_Treg_prop=table(Liver_FOXP3_Treg$orig.ident)/table(normal_liver$orig.ident)


Treg_data=data.frame(
  source=c(rep('CRLM',17),rep('HCC',8),rep('adj_normal',2+10)),  
  proportion=c(unname(CRLM_tumor.Treg.prop),unname(zzm_HCC_PT_Tcell.Treg.prop),c(unname(CRLM_adj.Treg.prop), unname(Liver_FOXP3_Treg_prop)))
)

Treg_data$source<-factor(Treg_data$source,levels = c('adj_normal','HCC','CRLM'), ordered=T) 


compaired=list(c('CRLM','adj_normal'),c('CRLM','HCC'))
ggplot(Treg_data,aes(x=source,y=proportion,color=source))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(colour = "black"),
        axis.text.x = element_text(angle = 30,vjust=0.5,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,colour = "black"))+
  labs(title = "Treg proportion comparison (LM)")+ylab('Treg proportion')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('ANT','HCC','CRLM'))+
  scale_colour_discrete(labels=c('ANT','HCC','CRLM'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("Fig S7C_4.pdf",width = 4,height = 6,units = "in")



########Fig S7C_3#######
Tcell.LM=Tcell %>% subset(lesions=='Liver_meta')
Tcell.LM.Th17.prop=table(Tcell.LM$cluster_name,Tcell.LM$orig.ident)['CD4_IL17A_Th17',]/table(Tcell.LM$orig.ident)

unname(Tcell.LM.Th17.prop) 

zzm_HCC_Th17=zzm_HCC_PT_Tcell %>% subset(IL17A> 0)
table(zzm_HCC_Th17$Sample)
table(zzm_HCC_PT_Tcell$Sample)

Tcell_LM_adj_Th17=Tcell_LM %>% subset(IL17A> 0)

normal_liver_Th17=normal_liver %>% subset(IL17A> 0)

Th17_adj_data=data.frame(
  source=c(rep('CRLM',17),rep('HCC',8),rep('adj_normal',12)),
  proportion=c(unname(Tcell.LM.Th17.prop),c(2/2747,3/1727,2/1624,5/856,0,0,0,0),c(10/4420,0,6/1818,10/6563,8/3010,20/4759,5/5014,1/1325,3/1705,0,0,0))
)


Th17_adj_data$source<-factor(Th17_adj_data$source,levels = c('adj_normal','HCC','CRLM'), ordered=T) 

compaired=list(c('CRLM','adj_normal'),c('CRLM','HCC'))

ggplot(Th17_adj_data,aes(x=source,y=proportion,color=source))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(colour = "black"),
        axis.text.x = element_text(angle = 30,vjust=0.5,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,colour = "black"))+
  labs(title = "Th17 proportion comparison (LM)")+ylab('Th17 proportion')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('ANT','HCC','CRLM'))+
  scale_colour_discrete(labels=c('ANT','HCC','CRLM'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("Fig S7C_3.pdf",width = 4,height = 6,units = "in")







