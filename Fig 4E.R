library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggsci)
library(ggpubr)
library(Seurat)
library(SeuratObject)

Tcell=readRDS("Tcell_in_ann_new_1224.rds")
Tcell=Tcell %>% subset(lesions=='Colon_tumor')
Meta_Th17=Tcell %>% subset(cluster_name=='CD4_IL17A_Th17')
Meta_Th17_TWEAK=Meta_Th17 %>% subset(TNFSF12>0)
Regev.CRC.intratumor.immune=readRDS("Regev_CRC_intratumor_immune.rds")
Regev.CRC.intratumor.Tcell=Regev.CRC.intratumor.immune %>% subset(clTopLevel==c('TNKILC'))
Nonmeta_Th17=Regev.CRC.intratumor.Tcell %>% subset(cl295v11SubFull=='cTNI05 (CD4+ IL17+)')
Nonmeta_Th17_TWEAK=Nonmeta_Th17 %>% subset(TNFSF12>0)

#####fig 4E right#### 
Meta_Th17_TWEAK$orig.ident=factor(Meta_Th17_TWEAK$orig.ident,levels=unique(Meta_Th17$orig.ident))
Meta_Th17$orig.ident=factor(Meta_Th17$orig.ident,levels=unique(Meta_Th17$orig.ident))

table(Meta_Th17_TWEAK$orig.ident)
table(Meta_Th17$orig.ident)

table(Meta_Th17_TWEAK$orig.ident)/table(Meta_Th17$orig.ident)

Meta_Th17_TWEAK_prop=as.data.frame(table(Meta_Th17_TWEAK$orig.ident)/table(Meta_Th17$orig.ident))
Meta_Th17_TWEAK_prop$Group='Meta'

Nonmeta_Th17_TWEAK$PID=factor(Nonmeta_Th17_TWEAK$PID,levels=unique(Nonmeta_Th17$PID))
Nonmeta_Th17$PID=factor(Nonmeta_Th17$PID,levels=unique(Nonmeta_Th17$PID))

table(Nonmeta_Th17_TWEAK$PID)
table(Nonmeta_Th17$PID)

table(Nonmeta_Th17_TWEAK$PID)/table(Nonmeta_Th17$PID)

Nonmeta_Th17_TWEAK_prop=as.data.frame(table(Nonmeta_Th17_TWEAK$PID)/table(Nonmeta_Th17$PID))
Nonmeta_Th17_TWEAK_prop$Group='Nonmeta'

Th17_TWEAK_prop=rbind(Meta_Th17_TWEAK_prop,Nonmeta_Th17_TWEAK_prop)

Th17_TWEAK_prop$Group=factor(Th17_TWEAK_prop$Group,levels = c('Nonmeta','Meta'),ordered = T)

compaired=list(c('Meta','Nonmeta'))
ggplot(Th17_TWEAK_prop,aes(x=Group,y=Freq,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+ylab('TNFSF12+ Th17 proprotion')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("fig 4E right.pdf",width = 3,height = 4.5,units = "in")




#####fig 4E left#### 
Meta_sub=Meta_Th17[rownames(Meta_Th17) %in% 'TNFSF12']

Idents(Meta_sub) = Meta_sub$orig.ident 
Meta_exp_by_sample = AverageExpression(Meta_sub, assays='RNA', slot='data')[[1]]
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample<-t(Meta_exp_by_sample)
head(Meta_exp_by_sample) 
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample$orig.ident<-rownames(Meta_exp_by_sample)
rownames(Meta_exp_by_sample)<-NULL
Meta_exp_by_sample$Group='Meta'

Non_Meta_sub=Nonmeta_Th17[rownames(Nonmeta_Th17) %in% 'TNFSF12']

Idents(Non_Meta_sub) = Non_Meta_sub$PID
Non_Meta_exp_by_sample = AverageExpression(Non_Meta_sub, assays='RNA', slot='data')[[1]]
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample<-t(Non_Meta_exp_by_sample)
head(Non_Meta_exp_by_sample) 
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample$orig.ident<-rownames(Non_Meta_exp_by_sample)
rownames(Non_Meta_exp_by_sample)<-NULL
Non_Meta_exp_by_sample$Group='Non_Meta'

df_exp=rbind(Meta_exp_by_sample,Non_Meta_exp_by_sample)

df_exp$Group=factor(df_exp$Group,levels = c('Non_Meta','Meta'),ordered = T)

compaired=list(c('Meta','Non_Meta'))
ggplot(df_exp,aes(x=Group,y=V1,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+ylab('TWEAK expression level in Th17 cells')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("fig 4E left.pdf",width = 3,height = 4.5,units = "in")





###LTA for revise####
Meta_sub=Meta_Th17[rownames(Meta_Th17) %in% 'LTA']

Idents(Meta_sub) = Meta_sub$orig.ident #
Meta_exp_by_sample = AverageExpression(Meta_sub, assays='RNA', slot='data')[[1]]
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample<-t(Meta_exp_by_sample)
head(Meta_exp_by_sample) 
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample$orig.ident<-rownames(Meta_exp_by_sample)
rownames(Meta_exp_by_sample)<-NULL
Meta_exp_by_sample$Group='Meta'

Non_Meta_sub=Nonmeta_Th17[rownames(Nonmeta_Th17) %in% 'LTA']

Idents(Non_Meta_sub) = Non_Meta_sub$PID
Non_Meta_exp_by_sample = AverageExpression(Non_Meta_sub, assays='RNA', slot='data')[[1]]
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample<-t(Non_Meta_exp_by_sample)
head(Non_Meta_exp_by_sample) 
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample$orig.ident<-rownames(Non_Meta_exp_by_sample)
rownames(Non_Meta_exp_by_sample)<-NULL
Non_Meta_exp_by_sample$Group='Non_Meta'

df_exp=rbind(Meta_exp_by_sample,Non_Meta_exp_by_sample)

df_exp$Group=factor(df_exp$Group,levels = c('Non_Meta','Meta'),ordered = T)

compaired=list(c('Meta','Non_Meta'))
ggplot(df_exp,aes(x=Group,y=V1,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+ylab('LTA expression level in Th17 cells')+xlab('')+
  #ggpubr::stat_compare_means(method = 'wilcox.test',alternative='greater')
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("LTA Th17 expression.pdf",width = 3,height = 4.5,units = "in")




#Treg
Meta_Treg=Tcell %>% subset(cluster_name=='CD4_FOXP3_Treg')
Nonmeta_Treg=Regev.CRC.intratumor.Tcell %>% subset(cl295v11SubFull=='cTNI08 (CD4+ Treg)')

Meta_sub=Meta_Treg[rownames(Meta_Treg) %in% 'LTA']

Idents(Meta_sub) = Meta_sub$orig.ident #
Meta_exp_by_sample = AverageExpression(Meta_sub, assays='RNA', slot='data')[[1]]
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample<-t(Meta_exp_by_sample)
head(Meta_exp_by_sample) 
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample$orig.ident<-rownames(Meta_exp_by_sample)
rownames(Meta_exp_by_sample)<-NULL
Meta_exp_by_sample$Group='Meta'

Non_Meta_sub=Nonmeta_Treg[rownames(Nonmeta_Treg) %in% 'LTA']

Idents(Non_Meta_sub) = Non_Meta_sub$PID
Non_Meta_exp_by_sample = AverageExpression(Non_Meta_sub, assays='RNA', slot='data')[[1]]
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample<-t(Non_Meta_exp_by_sample)
head(Non_Meta_exp_by_sample) 
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample$orig.ident<-rownames(Non_Meta_exp_by_sample)
rownames(Non_Meta_exp_by_sample)<-NULL
Non_Meta_exp_by_sample$Group='Non_Meta'

df_exp=rbind(Meta_exp_by_sample,Non_Meta_exp_by_sample)

df_exp$Group=factor(df_exp$Group,levels = c('Non_Meta','Meta'),ordered = T)

compaired=list(c('Meta','Non_Meta'))
ggplot(df_exp,aes(x=Group,y=V1,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+ylab('LTA expression level in Treg cells')+xlab('')+
  #ggpubr::stat_compare_means(method = 'wilcox.test',alternative='greater')
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
#ggsave("LTA Treg expression.pdf",width = 3,height = 4.5,units = "in")




