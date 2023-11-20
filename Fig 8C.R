library(Seurat)
library(dplyr)
library(ggsci)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(GSVA)


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


Regev.CRC.intratumor.monomacro=readRDS('Regev.CRC.intratumor.monomacro_ann.rds')
mono_macro=readRDS("myeloid_major_20220104_final.rds")
mono_macro.PT=mono_macro %>% subset(lesions=='Colon_tumor') 
mono_macro.LM=mono_macro %>% subset(lesions=='Liver_meta') 


Macro_c7_CXCL10=mono_macro.PT %>% subset(cluster_name=='Macro_c7_CXCL10')  
Macro_c7_CXCL10$PID_new=paste(Macro_c7_CXCL10$orig.ident, '_','Macro_c7_CXCL10')

Macro_c10_MT1G=mono_macro.PT %>% subset(cluster_name=='Macro_c10_MT1G')  
Macro_c10_MT1G$PID_new=paste(Macro_c10_MT1G$orig.ident, '_','Macro_c10_MT1G')

Macro_c1_CCL18=mono_macro.PT %>% subset(cluster_name=='Macro_c1_CCL18')  
Macro_c1_CCL18$PID_new=paste(Macro_c1_CCL18$orig.ident, '_','Macro_c1_CCL18')

Macro_c3_CD163L1=mono_macro.PT %>% subset(cluster_name=='Macro_c3_CD163L1')  
Macro_c3_CD163L1$PID_new=paste(Macro_c3_CD163L1$orig.ident, '_','Macro_c3_CD163L1')

Macro_c11_HSPH1=mono_macro.PT %>% subset(cluster_name=='Macro_c11_HSPH1')  
Macro_c11_HSPH1$PID_new=paste(Macro_c11_HSPH1$orig.ident, '_','Macro_c11_HSPH1')

Macro_c12_CD5L=mono_macro.PT %>% subset(cluster_name=='Macro_c12_CD5L')  
Macro_c12_CD5L$PID_new=paste(Macro_c12_CD5L$orig.ident, '_','Macro_c12_CD5L')

met_Macro_c3_CD163L1_LM=mono_macro.LM %>% subset(cluster_name=='Macro_c3_CD163L1')  
met_Macro_c3_CD163L1_LM$PID_new=paste(met_Macro_c3_CD163L1_LM$orig.ident, '_','met_Macro_c3_CD163L1_LM')

met_Macro_c1_CCL18_LM=mono_macro.LM %>% subset(cluster_name=='Macro_c1_CCL18')  
met_Macro_c1_CCL18_LM$PID_new=paste(met_Macro_c1_CCL18_LM$orig.ident, '_','met_Macro_c1_CCL18_LM')

met_Macro_c7_CXCL10_LM=mono_macro.LM %>% subset(cluster_name=='Macro_c7_CXCL10')  
met_Macro_c7_CXCL10_LM$PID_new=paste(met_Macro_c7_CXCL10_LM$orig.ident, '_','met_Macro_c7_CXCL10_LM')

met_Macro_c10_MT1G_LM=mono_macro.LM %>% subset(cluster_name=='Macro_c10_MT1G')  
met_Macro_c10_MT1G_LM$PID_new=paste(met_Macro_c10_MT1G_LM$orig.ident, '_','met_Macro_c10_MT1G_LM')

met_Macro_c11_HSPH1_LM=mono_macro.LM %>% subset(cluster_name=='Macro_c11_HSPH1')  
met_Macro_c11_HSPH1_LM$PID_new=paste(met_Macro_c11_HSPH1_LM$orig.ident, '_','met_Macro_c11_HSPH1_LM')

met_Macro_c12_CD5L_LM=mono_macro.LM %>% subset(cluster_name=='Macro_c12_CD5L')  
met_Macro_c12_CD5L_LM$PID_new=paste(met_Macro_c12_CD5L_LM$orig.ident, '_','met_Macro_c12_CD5L_LM')


cM12_Macro_MT1G=Regev.CRC.intratumor.monomacro %>% subset(cluster_name=='cM12_Macro_MT1G')  
cM12_Macro_MT1G$PID_new=paste(cM12_Macro_MT1G$PID, '_','cM12_Macro_MT1G')

cM11_Macro_CXCL10=Regev.CRC.intratumor.monomacro %>% subset(cluster_name=='cM11_Macro_CXCL10')  
cM11_Macro_CXCL10$PID_new=paste(cM11_Macro_CXCL10$PID, '_','cM11_Macro_CXCL10')

cM01_Macro_APOE=Regev.CRC.intratumor.monomacro %>% subset(cluster_name=='cM01_Macro_APOE')  
cM01_Macro_APOE$PID_new=paste(cM01_Macro_APOE$PID, '_','cM01_Macro_APOE')

cM02_Macro_SPP1=Regev.CRC.intratumor.monomacro %>% subset(cluster_name=='cM02_Macro_SPP1')  
cM02_Macro_SPP1$PID_new=paste(cM02_Macro_SPP1$PID, '_','cM02_Macro_SPP1')

cM03_Macro_C1QA=Regev.CRC.intratumor.monomacro %>% subset(cluster_name=='cM03_Macro_C1QA')  
cM03_Macro_C1QA$PID_new=paste(cM03_Macro_C1QA$PID, '_','cM03_Macro_C1QA')


Macro_merge=Macro_c7_CXCL10 %>% merge(Macro_c10_MT1G) %>% merge(Macro_c1_CCL18) %>% merge(Macro_c3_CD163L1) %>% merge(Macro_c11_HSPH1) %>% merge(Macro_c12_CD5L) %>%
  merge(met_Macro_c3_CD163L1_LM) %>% merge(met_Macro_c1_CCL18_LM) %>% merge(met_Macro_c7_CXCL10_LM) %>% 
  merge(met_Macro_c10_MT1G_LM) %>% merge(met_Macro_c11_HSPH1_LM) %>% merge(met_Macro_c12_CD5L_LM) %>% 
  merge(cM12_Macro_MT1G) %>% merge(cM11_Macro_CXCL10) %>% merge(cM01_Macro_APOE) %>% merge(cM02_Macro_SPP1) %>% merge(cM03_Macro_C1QA)


samplelist=unique(Macro_merge$PID_new)
Macro_merge.score=EMTscore(Macro_merge, samplelist)



Macro_EMTscore_data=data.frame(
  source=c(rep('met_Macro_c7_CXCL10_PT',length(unique(Macro_c7_CXCL10$PID_new))), #13
           rep('met_Macro_c10_MT1G_PT',length(unique(Macro_c10_MT1G$PID_new))),   #9
           rep('met_Macro_c1_CCL18_PT',length(unique(Macro_c1_CCL18$PID_new))),   #15
           rep('met_Macro_c3_CD163L1_PT',length(unique(Macro_c3_CD163L1$PID_new))), #16
           rep('met_Macro_c11_HSPH1_PT',length(unique(Macro_c11_HSPH1$PID_new))),   #13
           rep('met_Macro_c12_CD5L_PT',length(unique(Macro_c12_CD5L$PID_new))),   #4
           
           rep('met_Macro_c3_CD163L1_LM',length(unique(met_Macro_c3_CD163L1_LM$PID_new))), ###16
           rep('met_Macro_c1_CCL18_LM',length(unique(met_Macro_c1_CCL18_LM$PID_new))), #17
           rep('met_Macro_c7_CXCL10_LM',length(unique(met_Macro_c7_CXCL10_LM$PID_new))),  #16
           
           rep('met_Macro_c10_MT1G_LM',length(unique(met_Macro_c10_MT1G_LM$PID_new))),  #14
           rep('met_Macro_c11_HSPH1_LM',length(unique(met_Macro_c11_HSPH1_LM$PID_new))),  #13
           rep('met_Macro_c12_CD5L_LM',length(unique(met_Macro_c12_CD5L_LM$PID_new))),   #16
           
           rep('nonmet_cM12_Macro_MT1G',length(unique(cM12_Macro_MT1G$PID_new))), #26
           rep('nonmet_cM11_Macro_CXCL10',length(unique(cM11_Macro_CXCL10$PID_new))), #27
           rep('nonmet_cM01_Macro_APOE',length(unique(cM01_Macro_APOE$PID_new))),  #27
           rep('nonmet_cM02_Macro_SPP1',length(unique(cM02_Macro_SPP1$PID_new))),  #27
           rep('nonmet_cM03_Macro_C1QA',length(unique(cM03_Macro_C1QA$PID_new)))), #27
  EMTscore=c(unname(Macro_merge.score[,1:13]),
             unname(Macro_merge.score[,14:22]),
             unname(Macro_merge.score[,23:37]),
             unname(Macro_merge.score[,38:53]),
             unname(Macro_merge.score[,54:66]),
             unname(Macro_merge.score[,67:70]),#
             
             unname(Macro_merge.score[,71:86]),
             unname(Macro_merge.score[,87:103]),
             unname(Macro_merge.score[,104:119]),#+49
             
             unname(Macro_merge.score[,120:133]),
             unname(Macro_merge.score[,134:146]),
             unname(Macro_merge.score[,147:162]), 
             
             unname(Macro_merge.score[,163:188]),
             unname(Macro_merge.score[,189:215]),
             unname(Macro_merge.score[,216:242]),
             unname(Macro_merge.score[,243:269]),
             unname(Macro_merge.score[,270:296]))
)


plot.df1=Macro_EMTscore_data %>% subset(source%in%c("met_Macro_c3_CD163L1_PT","met_Macro_c1_CCL18_PT","met_Macro_c7_CXCL10_PT",
                                                    "met_Macro_c11_HSPH1_PT","met_Macro_c12_CD5L_PT","met_Macro_c10_MT1G_PT"))
plot.df1$source<-factor(plot.df1$source,levels = c("met_Macro_c3_CD163L1_PT","met_Macro_c1_CCL18_PT","met_Macro_c7_CXCL10_PT",
                                                   "met_Macro_c11_HSPH1_PT","met_Macro_c12_CD5L_PT","met_Macro_c10_MT1G_PT"), ordered=T) 


###Fig 8C####
ggplot(plot.df1,aes(x=source,y=EMTscore,color=source))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,vjust=0.9,hjust=0.9),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "Macrophage EMT relatedness (PT)")+ylab('EMT relatedness')+
  ggsignif::geom_signif(xmin=2,xmax=5, 
            y_position=0.28,  tip_length = 0.04,extend_line = 0,
            annotations = 'Wilcox p=3.686e-10',color='black')+
  geom_segment(aes(x=1,y=0.25,xend=3,yend=0.25),inherit.aes = F,color='black')+
  geom_segment(aes(x=4,y=0.25,xend=6,yend=0.25),inherit.aes = F,color='black')+xlab('')+
  scale_x_discrete(labels=c("Macro_c3_CD163L1","Macro_c1_CCL18","Macro_c7_CXCL10",
                            "Macro_c11_HSPH1","Macro_c12_CD5L","Macro_c10_MT1G"))+
  scale_colour_discrete(labels=c("Macro_c3_CD163L1","Macro_c1_CCL18","Macro_c7_CXCL10",
                                 "Macro_c11_HSPH1","Macro_c12_CD5L","Macro_c10_MT1G"),type = c(pal_d3()(6)))

ggsave("Fig 8C.pdf",width = 6,height = 5,units = "in")




