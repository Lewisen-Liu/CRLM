library(Seurat)
library(ggsci)
library(dplyr)
library(ggplot2)

####Fig 8A####
Myeloid.major=readRDS("myeloid_major_20220104_final.rds")

Myeloid.major$cluster_name=as.character(Myeloid.major$cluster_name)
Myeloid.major@meta.data$cluster_name[Myeloid.major@meta.data$cluster_name%in%c('Macro_c5_ANGPTL4')]='Mono_c5_ANGPTL4'

Myeloid.major$cluster_name <- factor(Myeloid.major$cluster_name, levels=c("Macro_c1_CCL18","Mono_c2_FCN1","Macro_c3_CD163L1","cDC2_c4_CLEC10A","Mono_c5_ANGPTL4","Cycling_Myeloid_c6",
                                                                          "Macro_c7_CXCL10","DC_c8_LAMP3","pDC_c9_LILRA4","Macro_c10_MT1G","Macro_c11_HSPH1","Macro_c12_CD5L","cDC1_c13_CLEC9A"), ordered=T)

mycolor=c('#46A5DA','#FAEC95','#5080CC','#EB4C20','#DEA077','#BF76CC','#9FCCA2','#40BFCC','#CC6F62','#8DADC4',
          '#CC3E74','#7373EA','#FEDDAA','#6BCC50','#BBB8CC','#D7DE5A','#8A52CC','#CC5323','#E288D1','#48CCAD',
          '#ACCCA3','#C936C2','#E55F25','#B9CEDB','#8F8FCC','#C19FCC','#F22727','#3975AC','#F5BC1D','#F1A9BA')
DimPlot(Myeloid.major,seed = 1, reduction = "umap",label = F,group.by = 'cluster_name',label.size = 6) +
  scale_color_manual(values = mycolor)
ggsave("Fig 8A.pdf",width = 5.5,height = 4,units = "in")


####Fig S9B####
DimPlot(Myeloid.major,seed = 1, reduction = "umap",label = F,group.by = 'lesions',label.size = 6) +scale_color_igv()
ggsave("Fig S9B.pdf",width = 5,height = 4,units = "in")



####Fig S9A####
FeaturePlot(Myeloid.major, features=c('CCL18','FCN1','CD163L1',
                                      'CLEC10A','ANGPTL4','MKI67',
                                      'CXCL10','LAMP3','LILRA4','MT1G','HSPH1','CD5L','CLEC9A'),reduction = "umap",cols = c("#DFDFDF","#FF0000"),ncol = 5) 
ggsave("Fig S9A.pdf",width = 12,height = 7,units = "in")



