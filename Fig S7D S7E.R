library(ggplot2)
library(ggsignif)
library(dplyr)
library(ggsci)
library(ggpubr)
library(Seurat)
library(SeuratObject)

Tcell=readRDS("Tcell_in_ann_new_1224.rds")
#Tcell=Tcell %>% subset(lesions=='Colon_tumor')
Meta_Th17=Tcell %>% subset(cluster_name=='CD4_IL17A_Th17')
#Meta_Th17_TWEAK=Meta_Th17 %>% subset(TNFSF12>0)
Regev.CRC.intratumor.immune=readRDS("Regev_CRC_intratumor_immune.rds")
Regev.CRC.intratumor.Tcell=Regev.CRC.intratumor.immune %>% subset(clTopLevel==c('TNKILC'))
Nonmeta_Th17=Regev.CRC.intratumor.Tcell %>% subset(cl295v11SubFull=='cTNI05 (CD4+ IL17+)')
#Nonmeta_Th17_TWEAK=Nonmeta_Th17 %>% subset(TNFSF12>0)


#####Fig S7D Th17 score: mCRC vs nmCRC####

Meta_Th17$Group='mCRC'
Nonmeta_Th17$Group='nmCRC'
merge.Th17=Meta_Th17 %>% merge(Nonmeta_Th17)

Th17.list<-list(c('IL17A','CCL20','RORA','STAT3','IL21','IL23R','RORC'))
merge.Th17=NormalizeData(merge.Th17,normalization.method = "LogNormalize",scale.factor = 10000)

merge.Th17=ScaleData(merge.Th17,features = c('IL17A','CCL20','RORA','STAT3','IL21','IL23R','RORC'))
merge.Th17<-AddModuleScore(object = merge.Th17,features = Th17.list, name='Th17') 
colnames(merge.Th17@meta.data)[colnames(merge.Th17@meta.data)=='Th171']='Th17'

meta.Th17.df=merge.Th17 %>% subset(Group%in%'mCRC' & lesions=='Colon_tumor')
Meta_exp_by_sample=meta.Th17.df@meta.data %>% group_by(orig.ident) %>% summarise(Th17score=mean(Th17, na.rm = TRUE))
Meta_exp_by_sample$Group='Meta'
nonmeta.Th17.df=merge.Th17 %>% subset(Group%in%'nmCRC')
Non_Meta_exp_by_sample=nonmeta.Th17.df@meta.data %>% group_by(orig.ident) %>% summarise(Th17score=mean(Th17, na.rm = TRUE))
Non_Meta_exp_by_sample$Group='Non_Meta'

df_exp=rbind(Meta_exp_by_sample,Non_Meta_exp_by_sample)

df_exp$Group=factor(df_exp$Group,levels = c('Non_Meta','Meta'),ordered = T)

compaired=list(c('Meta','Non_Meta'))
ggplot(df_exp,aes(x=Group,y=Th17score,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(position=position_jitter(width = 0.1, height=0),alpha=0.5,size=2)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  labs(title = "")+ylab('Th17 score')+xlab('')+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.1, map_signif_level = F, test = 'wilcox.test',color='black')+
  scale_x_discrete(labels=c('nmCRC','mCRC'))+
  scale_colour_discrete(labels=c('nmCRC','mCRC'),type = c('#1f77b4','#ff7f0e','#2ca02c'))

ggsave("Fig S7D Th17score_mCRC_vs_nmCRC.pdf",width = 3,height = 4.5,units = "in")








#####Fig S7E Th17 score: PT vs LM####
Th17_score<-merge.Th17@meta.data %>% group_by(orig.ident) %>% dplyr::summarise(Th17_score=mean(Th17))
head(merge.Th17@meta.data)
clin<-merge.Th17@meta.data[,c('orig.ident','patient','lesions')] %>% unique()
clin
rownames(clin)<-NULL

Th17_score<-left_join(Th17_score,clin,by='orig.ident')
Th17_score=Th17_score%>%arrange(desc(lesions))

Th17_score$group=c(rep('mCRC_LM',17),rep('mCRC_PT',16),rep('nmCRC',26))
Th17_score$group=factor(Th17_score$group,levels = c('nmCRC','mCRC_PT','mCRC_LM'),ordered = T)

plot.data<-Th17_score[1:33,]
head(plot.data)
plot.data.pared<-plot.data[-c(3,5,13,16,20,31,33),] 
plot.data.pared<-plot.data.pared %>% arrange(lesions,patient)  
plot.data.pared<-plot.data.pared[-c(9,22),] 
wilcox.test(plot.data.pared$Th17_score[1:12],plot.data.pared$Th17_score[13:24],paired = T,alternative = c("greater")) 

ggpaired(plot.data.pared, x = 'lesions', y = 'Th17_score',
         color='lesions',fill = 'lesions', line.color = "#7f7f7f",linetype ='dashed', line.size = 0.4, point.size = 2)+   
  scale_color_igv()+
  scale_fill_manual(values = c('#5050FF80','#CE3D3280'))+
  labs(x = '', y = 'Th17_score', title = '', subtitle = '')+geom_signif(xmin=1,
                                                                       xmax=2, 
                                                                       y_position=0.38,
                                                                       annotations = 'p=0.0034')

ggsave("Fig S7E Th17score_PT_vs_LM.pdf",width = 2.5,height = 4,units = "in")



