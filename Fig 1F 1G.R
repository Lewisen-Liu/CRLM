#figure 1F, 1G
library(immunedeconv)
library(dplyr)
library(stringi)
library(biomaRt) #ensembl
library(stringr)
library(xCell)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library(GSVA)

tpms=readRDS('bulk_data_Fig1G_tpms.rds')
merge_df=readRDS("bulk_data_Fig1G_counts.rds")

Th17=as.data.frame(tpms[c('IL17A','IL21','IL23R','RORA','STAT3','CCL20','RORC'),]) 
Th17=as.data.frame(t(Th17))

Th17$group='mCRC'
Th17$group[34:226]='nmCRC'
Th17$group=factor(Th17$group,levels = c('nmCRC','mCRC'),ordered = T)

compaired=list(c('mCRC','nmCRC'))

#Fig 1F#######
p1=ggplot(Th17,aes(x=group,y=IL21,color=group))+
  geom_violin()+
  scale_color_igv()+
  geom_jitter(position=position_jitter(width = 0.2, height=0),alpha=1,size=0.5)+  
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0, hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+xlab('')+ylab('Expression level of IL21 (TPM)')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')

p2=ggplot(Th17,aes(x=group,y=IL23R,color=group))+
  geom_violin()+
  scale_color_igv()+
  geom_jitter(position=position_jitter(width = 0.2, height=0),alpha=1,size=0.5)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0, hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+xlab('')+ylab('Expression level of IL23R (TPM)')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')

p3=ggplot(Th17,aes(x=group,y=RORA,color=group))+
  geom_violin()+
  scale_color_igv()+
  geom_jitter(position=position_jitter(width = 0.2, height=0),alpha=1,size=0.5)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0, hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+xlab('')+ylab('Expression level of RORA (TPM)')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')

p4=ggplot(Th17,aes(x=group,y=STAT3,color=group))+
  geom_violin()+
  scale_color_igv()+
  geom_jitter(position=position_jitter(width = 0.2, height=0),alpha=1,size=0.5)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0, hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+xlab('')+ylab('Expression level of STAT3 (TPM)')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')

p3+p1+p4+p2+patchwork::plot_layout(ncol = 4, guides = 'collect')
ggsave("Fig 1F.pdf",width = 8,height = 5,units = "in")
ggsave("Fig 1F.pdf",width = 8,height = 3,units = "in")



#Fig 1G######
dim(merge_df)
merge_df[1:4,1:4]

res.matrix = gsva(as.matrix(merge_df[,1:226]), list(c('IL17A','CCL20','RORA','STAT3','IL21','IL23R','RORC')), method='ssgsea', kcdf="Poisson",parallel.sz=10)

res.matrix[1:4]
wilcox.test(res.matrix[1:33],res.matrix[34:226])

Th17=as.data.frame(t(res.matrix))
colnames(Th17)='Th17'

wilcox.test(Th17$Th17[1:33],Th17$Th17[34:226])
mean(Th17$Th17[1:33])
mean(Th17$Th17[34:226])

Th17$group='mCRC'
Th17$group[34:226]='nmCRC'
Th17$group=factor(Th17$group,levels = c('nmCRC','mCRC'),ordered = T)

compaired=list(c('mCRC','nmCRC'))
ggplot(Th17,aes(x=group,y=Th17,color=group))+
  geom_violin()+
  scale_color_igv()+
  geom_jitter(position=position_jitter(width = 0.3, height=0),alpha=1,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0, hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+xlab('')+ylab('Th17 score (xCell)')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')

ggsave("Fig 1G.pdf",width = 4,height = 5,units = "in")

