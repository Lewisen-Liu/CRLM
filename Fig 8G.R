library(Seurat)
library(ggsci)

Regev.CRC.intratumor.monomacro=readRDS('Regev.CRC.intratumor.monomacro_ann.rds')
mono_macro=readRDS("myeloid_major_20220104_final.rds")
mono_macro.PT=mono_macro %>% subset(lesions=='Colon_tumor') 
mono_macro.LM=mono_macro %>% subset(lesions=='Liver_meta') 

Meta=mono_macro.PT
Non_Meta=Regev.CRC.intratumor.monomacro

Meta_sub=Meta[rownames(Meta) %in% 'CD163L1']
Idents(Meta_sub) = Meta_sub$orig.ident #CD163L1
Meta_exp_by_sample = AverageExpression(Meta_sub, assays='RNA', slot='data')[[1]]
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample<-t(Meta_exp_by_sample)
head(Meta_exp_by_sample) 
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample$orig.ident<-rownames(Meta_exp_by_sample)
rownames(Meta_exp_by_sample)<-NULL
Meta_exp_by_sample$Group='Meta'

Non_Meta_sub=Non_Meta[rownames(Non_Meta) %in% 'CD163L1']
Idents(Non_Meta_sub) = Non_Meta_sub$PID
Non_Meta_exp_by_sample = AverageExpression(Non_Meta_sub, assays='RNA', slot='data')[[1]]
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample<-t(Non_Meta_exp_by_sample)
head(Non_Meta_exp_by_sample) 
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample$orig.ident<-rownames(Non_Meta_exp_by_sample)
rownames(Non_Meta_exp_by_sample)<-NULL
Non_Meta_exp_by_sample$Group='Non_Meta'

CD163L1_exp=rbind(Meta_exp_by_sample,Non_Meta_exp_by_sample)
CD163L1_exp$Group=factor(CD163L1_exp$Group,levels = c('Non_Meta','Meta'),ordered = T)

compaired=list(c('Meta','Non_Meta'))

##Fig 8G######
ggplot(CD163L1_exp,aes(x=Group,y=V1,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+ylab('CD163L1 expression level in Macrophage')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("Fig 8G.pdf",width = 3,height = 4,units = "in")



