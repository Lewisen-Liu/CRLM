library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(Seurat)
Regev.CRC.intratumor.monomacro=readRDS('Regev.CRC.intratumor.monomacro_ann.rds')
mono_macro=readRDS("myeloid_major_20220104_final.rds")

#####Fig 8J#####
Regev.macro=Regev.CRC.intratumor.monomacro %>% subset(cluster_name%in%c('cM01_Macro_APOE','cM03_Macro_C1QA','cM11_Macro_CXCL10','cM12_Macro_MT1G')) #'cM02_Macro_SPP1'
mCRC.macro=mono_macro %>% subset(cluster_name%in%c('Macro_c5_ANGPTL4','Macro_c3_CD163L1','Macro_c1_CCL18','Macro_c11_HSPH1','Macro_c7_CXCL10','Macro_c10_MT1G','Macro_c12_CD5L'))
Regev.macro$group='nmCRC'
mCRC.macro$group='mCRC'

merge.macro=merge(Regev.macro,mCRC.macro)
merge.macro$group=factor(merge.macro$group,levels = c('nmCRC','mCRC'),ordered = T)

VlnPlot(merge.macro,features = c('SPP1'),group.by = 'group',pt.size = 0)+ggpubr::stat_compare_means(method = "t.test",label.x = 1.2)+scale_fill_igv()+xlab('')+
  theme(axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5))
ggsave("Fig 8J.pdf",width = 4,height = 4,units = "in")




####Fig 8K#####
mCRC_SPP1=mono_macro %>% subset(SPP1>0)
nmCRC_SPP1=Regev.CRC.intratumor.monomacro %>% subset(cluster_name%in%c('cM02_Macro_SPP1'))

dim(mCRC_SPP1)[2]/dim(mono_macro)[2]
dim(nmCRC_SPP1)[2]/dim(Regev.CRC.intratumor.monomacro)[2]

mCRC_SPP1$orig.ident=factor(mCRC_SPP1$orig.ident,levels=unique(mono_macro$orig.ident))
mono_macro$orig.ident=factor(mono_macro$orig.ident,levels=unique(mono_macro$orig.ident))

table(mCRC_SPP1$orig.ident)
table(mono_macro$orig.ident)

table(mCRC_SPP1$orig.ident)/table(mono_macro$orig.ident)

mCRC_SPP1_prop=as.data.frame(table(mCRC_SPP1$orig.ident)/table(mono_macro$orig.ident))
mCRC_SPP1_prop$Group='mono_macro'


nmCRC_SPP1$PID=factor(nmCRC_SPP1$PID,levels=unique(Regev.CRC.intratumor.monomacro$PID))
Regev.CRC.intratumor.monomacro$PID=factor(Regev.CRC.intratumor.monomacro$PID,levels=unique(Regev.CRC.intratumor.monomacro$PID))

table(nmCRC_SPP1$PID)
table(Regev.CRC.intratumor.monomacro$PID)

table(nmCRC_SPP1$PID)/table(Regev.CRC.intratumor.monomacro$PID)

nmCRC_SPP1_prop=as.data.frame(table(nmCRC_SPP1$PID)/table(Regev.CRC.intratumor.monomacro$PID))
nmCRC_SPP1_prop$Group='Regev.CRC.intratumor.monomacro'

SPP1_prop=rbind(mCRC_SPP1_prop,nmCRC_SPP1_prop)
SPP1_prop$Group=factor(SPP1_prop$Group,levels = c('Regev.CRC.intratumor.monomacro','mono_macro'),ordered = T)

compaired=list(c('mono_macro','Regev.CRC.intratumor.monomacro'))
ggplot(SPP1_prop,aes(x=Group,y=Freq,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  #scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5), #angle = 45,vjust=0.9,hjust=0.9
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+ylab('SPP1+ Macrophage ratio')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e'))
ggsave("Fig 8K.pdf",width = 3,height = 4,units = "in")






