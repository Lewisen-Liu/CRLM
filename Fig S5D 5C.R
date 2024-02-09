library(Seurat)
library(dplyr)
library(ggsci)
library(ggsignif)
library(ggplot2)
library(ggpubr)

######Fig S5D#######
Regev.CRC.intratumor.epi=readRDS("Regev_CRC_intratumor_tumor.rds")
Regev.CRC.intratumor.epi$group='nmCRC'

Tumor = readRDS("tumors_filter_by_cnv_man.rds")
Tumor.PT=Tumor %>% subset(lesions=='Colon_tumor')
Tumor.PT$group='mCRC_PT'
Tumor.LM=Tumor %>% subset(lesions=='Liver_meta')
Tumor.LM$group='mCRC_LM'

merge.tumor=Regev.CRC.intratumor.epi %>% merge(Tumor.PT) %>% merge(Tumor.LM)

IG_name<-rownames(merge.tumor)[grep("^IG[HLK]",rownames(merge.tumor))]
HSP_name<-rownames(merge.tumor)[grep("^HSP",rownames(merge.tumor))]
keep_name<-setdiff(rownames(merge.tumor), c(IG_name,"JCHAIN")) %>% setdiff(HSP_name)
merge.tumor<-subset(merge.tumor, features=keep_name)

excludeGene<-read.table("excludeGene.txt",header=TRUE)
merge.tumor<-merge.tumor[!rownames(merge.tumor) %in% excludeGene$Gene]

EMT.list<-list(c('TIMP1','LAMC2','TGFBI','LAMA3','SDC4','PMEPA1','SPP1','ECM1','IGFBP4','CD59','FLNA','LGALS1','SERPINE2','VIM','PLOD3'))
merge.tumor=NormalizeData(merge.tumor,normalization.method = "LogNormalize",scale.factor = 10000)
EMT<-read.delim('EMT.txt',header = T)
EMT=EMT[-1,]
merge.tumor=ScaleData(merge.tumor,features = EMT)
merge.tumor<-AddModuleScore(object = merge.tumor,features = EMT.list, name='EMT') 
colnames(merge.tumor@meta.data)[colnames(merge.tumor@meta.data)=='EMT1']='EMT'

merge.tumor$group=factor(merge.tumor$group,levels = c('nmCRC','mCRC_PT','mCRC_LM'),ordered = T)
compaired=list(c('nmCRC','mCRC_PT'))

colors <- colorRampPalette(c("pink", "darkred"))(100)
colors[c(1,17,100)]

VlnPlot(merge.tumor, features = 'EMT',pt.size = 0, group.by = "group")+
  scale_fill_manual(values = colors[c(1,17,100)])+
  geom_boxplot(width=0.2,fill='white',outlier.shape = NA)+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0, hjust=0.5))+
  labs(title = "EMT score (cell level)")+xlab('')+ylab('EMT score')

ggsave("Fig S5D.pdf",width = 5,height = 5,units = "in")



#####Fig 5C####
EMT_score<-merge.tumor@meta.data %>% group_by(orig.ident) %>% dplyr::summarise(EMT_score=mean(EMT))
head(Tumor@meta.data)
clin<-Tumor@meta.data[,c('orig.ident','patient','lesions')] %>% unique()
clin
rownames(clin)<-NULL

EMT_score<-left_join(EMT_score,clin,by='orig.ident')
EMT_score=EMT_score%>%arrange(desc(lesions))

EMT_score$group=c(rep('mCRC_LM',17),rep('mCRC_PT',16),rep('nmCRC',27))
EMT_score$group=factor(EMT_score$group,levels = c('nmCRC','mCRC_PT','mCRC_LM'),ordered = T)

plot.data<-EMT_score[1:33,]
head(plot.data)
plot.data.pared<-plot.data[-c(3,5,13,16,20,31,33),] 
plot.data.pared<-plot.data.pared %>% arrange(lesions,patient)  
plot.data.pared<-plot.data.pared[-c(9,22),] 
wilcox.test(plot.data.pared$EMT_score[1:12],plot.data.pared$EMT_score[13:24],paired = T,alternative = c("less")) 

ggpaired(plot.data.pared, x = 'lesions', y = 'EMT_score',
         color='lesions',fill = 'lesions', line.color = "#7f7f7f",linetype ='dashed', line.size = 0.4, point.size = 2)+   
  scale_color_igv()+
  scale_fill_manual(values = c('#5050FF80','#CE3D3280'))+
  labs(x = '', y = 'EMT_score', title = '', subtitle = '')+geom_signif(xmin=1,
                                                                       xmax=2, 
                                                                       y_position=0.32,
                                                                       annotations = 'p=0.0034')

ggsave("Fig 5C.pdf",width = 2.5,height = 4,units = "in")



