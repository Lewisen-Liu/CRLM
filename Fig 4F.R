library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(Seurat)

Meta = readRDS("tumors_filter_by_cnv_man.rds")
Meta=NormalizeData(Meta,normalization.method = "LogNormalize",scale.factor = 10000)

Non_Meta = readRDS("Regev_CRC_intratumor_tumor.rds")
Non_Meta=NormalizeData(Non_Meta,normalization.method = "LogNormalize",scale.factor = 10000)


#####Fig 4F left########
Meta_sub=Meta[rownames(Meta) %in% 'TNFRSF12A']

Idents(Meta_sub) = Meta_sub$orig.ident 
Meta_exp_by_sample = AverageExpression(Meta_sub, assays='RNA', slot='data')[[1]]
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample<-t(Meta_exp_by_sample)
head(Meta_exp_by_sample)
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample$orig.ident<-rownames(Meta_exp_by_sample)
rownames(Meta_exp_by_sample)<-NULL
Meta_exp_by_sample$Group='Meta'

Non_Meta_sub=Non_Meta[rownames(Non_Meta) %in% 'TNFRSF12A']

Idents(Non_Meta_sub) = Non_Meta_sub$PID
Non_Meta_exp_by_sample = AverageExpression(Non_Meta_sub, assays='RNA', slot='data')[[1]]
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample<-t(Non_Meta_exp_by_sample)
head(Non_Meta_exp_by_sample) 
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample$orig.ident<-rownames(Non_Meta_exp_by_sample)
rownames(Non_Meta_exp_by_sample)<-NULL
Non_Meta_exp_by_sample$Group='Non_Meta'


Fn14_exp=rbind(Meta_exp_by_sample,Non_Meta_exp_by_sample)
Fn14_exp$Group=factor(Fn14_exp$Group,levels = c('Non_Meta','Meta'),ordered = T)
compaired=list(c('Meta','Non_Meta'))
ggplot(Fn14_exp,aes(x=Group,y=V1,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+ylab('TNFRSF12 expression level in tumor cells')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("Fig 4F left.pdf",width = 3,height = 4.5,units = "in")


#####Fig 4F right########
Meta_sub=Meta[rownames(Meta) %in% 'TNFRSF25']

Idents(Meta_sub) = Meta_sub$orig.ident 
Meta_exp_by_sample = AverageExpression(Meta_sub, assays='RNA', slot='data')[[1]]
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample<-t(Meta_exp_by_sample)
head(Meta_exp_by_sample)
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample$orig.ident<-rownames(Meta_exp_by_sample)
rownames(Meta_exp_by_sample)<-NULL
Meta_exp_by_sample$Group='Meta'



Non_Meta_sub=Non_Meta[rownames(Non_Meta) %in% 'TNFRSF25']

Idents(Non_Meta_sub) = Non_Meta_sub$PID
Non_Meta_exp_by_sample = AverageExpression(Non_Meta_sub, assays='RNA', slot='data')[[1]]
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample<-t(Non_Meta_exp_by_sample)
head(Non_Meta_exp_by_sample) 
Non_Meta_exp_by_sample = as.data.frame(Non_Meta_exp_by_sample)
Non_Meta_exp_by_sample$orig.ident<-rownames(Non_Meta_exp_by_sample)
rownames(Non_Meta_exp_by_sample)<-NULL
Non_Meta_exp_by_sample$Group='Non_Meta'



Fn14_exp=rbind(Meta_exp_by_sample,Non_Meta_exp_by_sample)
Fn14_exp$Group=factor(Fn14_exp$Group,levels = c('Non_Meta','Meta'),ordered = T)
compaired=list(c('Meta','Non_Meta'))
ggplot(Fn14_exp,aes(x=Group,y=V1,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+ylab('TNFRSF25 expression level in tumor cells')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("Fig 4F left.pdf",width = 3,height = 4.5,units = "in")









#####LTBR for revise########
Meta_sub=Meta[rownames(Meta) %in% 'LTBR']

Idents(Meta_sub) = Meta_sub$orig.ident 
Meta_exp_by_sample = AverageExpression(Meta_sub, assays='RNA', slot='data')[[1]]
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample<-t(Meta_exp_by_sample)
head(Meta_exp_by_sample) 
Meta_exp_by_sample = as.data.frame(Meta_exp_by_sample)
Meta_exp_by_sample$orig.ident<-rownames(Meta_exp_by_sample)
rownames(Meta_exp_by_sample)<-NULL
Meta_exp_by_sample$Group='Meta'



Non_Meta_sub=Non_Meta[rownames(Non_Meta) %in% 'LTBR']

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
  labs(title = "")+ylab('LTBR expression level in tumor cells')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave("LTBR tumor cell expression.pdf",width = 3,height = 4.5,units = "in")






