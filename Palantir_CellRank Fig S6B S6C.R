library(Seurat)
library(dplyr)
library(ggsci)
library(ggplot2)
library(SeuratDisk)


EMT<-read.delim('EMT.txt')
EMT<-EMT[-1,]

Tumor = readRDS("tumors_filter_by_cnv_man.rds")

IG_name<-rownames(Tumor)[grep("^IG[HLK]",rownames(Tumor))]
HSP_name<-rownames(Tumor)[grep("^HSP",rownames(Tumor))]
keep_name<-setdiff(rownames(Tumor), c(IG_name,"JCHAIN")) %>% setdiff(HSP_name)
Tumor<-subset(Tumor, features=keep_name)

excludeGene<-read.table("excludeGene.txt",header=TRUE)
Tumor<-Tumor[!rownames(Tumor) %in% excludeGene$Gene]


set.seed(1)

EMT.list<-list(EMT)
Tumor<-AddModuleScore(object = Tumor,features = EMT.list, name='EMT')
VlnPlot(Tumor, features = 'EMT1',pt.size = 0, group.by = "lesions")+
  scale_fill_igv()+
  geom_boxplot(width=0.2,fill='white')+
  ggpubr::stat_compare_means(method = 't.test',hjust = -0.25, label="p.format")

Tumor$EMT_score<-Tumor$EMT1




Pt02<-Tumor %>% subset(patient%in%'Pt02')
SaveH5Seurat(Pt02, filename = "Pt02.h5seurat", overwrite=TRUE)
Convert("Pt02.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt01<-Tumor %>% subset(patient%in%'Pt01')
SaveH5Seurat(Pt01, filename = "Pt01.h5seurat", overwrite=TRUE)
Convert("Pt01.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt03<-Tumor %>% subset(patient%in%'Pt03')
SaveH5Seurat(Pt03, filename = "Pt03.h5seurat", overwrite=TRUE)
Convert("Pt03.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt04<-Tumor %>% subset(patient%in%'Pt04')
SaveH5Seurat(Pt04, filename = "Pt04.h5seurat", overwrite=TRUE)
Convert("Pt04.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt05<-Tumor %>% subset(patient%in%'Pt05')
SaveH5Seurat(Pt05, filename = "Pt05.h5seurat", overwrite=TRUE)
Convert("Pt05.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt06<-Tumor %>% subset(patient%in%'Pt06')
SaveH5Seurat(Pt06, filename = "Pt06.h5seurat", overwrite=TRUE)
Convert("Pt06.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt07<-Tumor %>% subset(patient%in%'Pt07')
SaveH5Seurat(Pt07, filename = "Pt07.h5seurat", overwrite=TRUE)
Convert("Pt07.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt08<-Tumor %>% subset(patient%in%'Pt08')
SaveH5Seurat(Pt08, filename = "Pt08.h5seurat", overwrite=TRUE)
Convert("Pt08.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt10<-Tumor %>% subset(patient%in%'Pt10')
SaveH5Seurat(Pt10, filename = "Pt10.h5seurat", overwrite=TRUE)
Convert("Pt10.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt11<-Tumor %>% subset(patient%in%'Pt11')
SaveH5Seurat(Pt11, filename = "Pt11.h5seurat", overwrite=TRUE)
Convert("Pt11.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt15<-Tumor %>% subset(patient%in%'Pt15')
SaveH5Seurat(Pt15, filename = "Pt15.h5seurat", overwrite=TRUE)
Convert("Pt15.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt17<-Tumor %>% subset(patient%in%'Pt17')
SaveH5Seurat(Pt17, filename = "Pt17.h5seurat", overwrite=TRUE)
Convert("Pt17.h5seurat", dest = "h5ad", overwrite = TRUE)

Pt18<-Tumor %>% subset(patient%in%'Pt18')
SaveH5Seurat(Pt18, filename = "Pt18.h5seurat", overwrite=TRUE)
Convert("Pt18.h5seurat", dest = "h5ad", overwrite = TRUE)

#Then use the h5ad files as input for 'palantir_run.ipynb' and 'CellRank.ipynb'




####Palantir######
Pt01_res=read.csv("Pt01_annotations.csv")
Pt02_res=read.csv("Pt02_annotations.csv")
Pt03_res=read.csv("Pt03_annotations.csv")
Pt04_res=read.csv("Pt04_annotations.csv")
Pt05_res=read.csv("Pt05_annotations.csv")
Pt06_res=read.csv("Pt06_annotations.csv")
Pt07_res=read.csv("Pt07_annotations.csv")
Pt08_res=read.csv("Pt08_annotations.csv")
Pt10_res=read.csv("Pt10_annotations.csv")
Pt11_res=read.csv("Pt11_annotations.csv")
Pt15_res=read.csv("Pt15_annotations.csv")
Pt17_res=read.csv("Pt17_annotations.csv")
Pt18_res=read.csv("Pt18_annotations.csv")




res=Pt01_res %>% rbind(Pt02_res) %>% rbind(Pt03_res) %>% rbind(Pt04_res) %>% rbind(Pt05_res) %>% 
  rbind(Pt06_res) %>% rbind(Pt07_res) %>% rbind(Pt08_res) %>% rbind(Pt10_res) %>% rbind(Pt11_res) %>% 
  rbind(Pt15_res) %>% rbind(Pt17_res) %>% rbind(Pt18_res)

head(res)
res_palantir=res %>% group_by(patient,lesions) %>% summarise(palantir=mean(palantir_pseudotime, na.rm = TRUE))
res_palantir<-res_palantir %>% arrange(lesions,patient) 

wilcox.test(res_palantir$palantir[1:13],res_palantir$palantir[14:26],paired = T,alternative = c("less")) 

ggpaired(res_palantir, x = 'lesions', y = 'palantir',
         color='lesions',fill = 'lesions', line.color = "#7f7f7f",linetype ='dashed', line.size = 0.4, point.size = 2)+   
  scale_color_igv()+
  scale_fill_manual(values = c('#5050FF80','#CE3D3280'))+
  labs(x = '', y = 'Palantir pseudotime', title = '', subtitle = '')+geom_signif(xmin=1,
                                                                        xmax=2, 
                                                                        y_position=1.25,
                                                                        annotations = 'p=0.00037')
ggsave("Palantir_pseudotime_PT_vs_LM.pdf",width = 2.5,height = 4,units = "in")



pData.new=as.data.frame(matrix(nrow=0, ncol = 6))
colnames(pData.new)=c(colnames(res),'Pseudotime')
pt.list=unique(res$patient)
for (i in pt.list) {
  pD.tmp=subset(res, patient %in% i)
  pD.tmp$Pseudotime=(pD.tmp$palantir_pseudotime-min(pD.tmp$palantir_pseudotime))/(max(pD.tmp$palantir_pseudotime)-min(pD.tmp$palantir_pseudotime))*100 
  pData.new=rbind(pData.new, pD.tmp) 
}

pData.new$phase[pData.new$Pseudotime<=100/6]='E'
pData.new$phase[pData.new$Pseudotime>100/6 & pData.new$Pseudotime<=200/6]='H1'
pData.new$phase[pData.new$Pseudotime>200/6 & pData.new$Pseudotime<=300/6]='H2'
pData.new$phase[pData.new$Pseudotime>300/6 & pData.new$Pseudotime<=400/6]='H3'
pData.new$phase[pData.new$Pseudotime>400/6 & pData.new$Pseudotime<=500/6]='H4'
pData.new$phase[pData.new$Pseudotime>500/6 & pData.new$Pseudotime<=100]='M'



EMT_compare_E=c()
EMT_compare_H1=c()
EMT_compare_H2=c()
EMT_compare_H3=c()
EMT_compare_H4=c()
EMT_compare_M=c()

for (i in pt.list) {
  pD.tmp=subset(pData.new, patient%in%i)
  pD.E =pD.tmp%>% arrange(Pseudotime) %>% subset(Pseudotime<unname(quantile(pD.tmp$Pseudotime,1/6))) %>% summarise(EMT_score_new=mean(EMT_score))
  pD.H1 =pD.tmp%>% arrange(Pseudotime) %>% subset(Pseudotime<unname(quantile(pD.tmp$Pseudotime,2/6))) %>% subset(Pseudotime>unname(quantile(pD.tmp$Pseudotime,1/6))) %>% summarise(EMT_score_new=mean(EMT_score))
  pD.H2 =pD.tmp%>% arrange(Pseudotime) %>% subset(Pseudotime<unname(quantile(pD.tmp$Pseudotime,3/6))) %>% subset(Pseudotime>unname(quantile(pD.tmp$Pseudotime,2/6))) %>% summarise(EMT_score_new=mean(EMT_score))
  pD.H3 =pD.tmp%>% arrange(Pseudotime) %>% subset(Pseudotime<unname(quantile(pD.tmp$Pseudotime,4/6))) %>% subset(Pseudotime>unname(quantile(pD.tmp$Pseudotime,3/6))) %>% summarise(EMT_score_new=mean(EMT_score))
  pD.H4 =pD.tmp%>% arrange(Pseudotime) %>% subset(Pseudotime<unname(quantile(pD.tmp$Pseudotime,5/6))) %>% subset(Pseudotime>unname(quantile(pD.tmp$Pseudotime,4/6))) %>% summarise(EMT_score_new=mean(EMT_score))
  pD.M =pD.tmp%>% arrange(Pseudotime) %>% subset(Pseudotime>unname(quantile(pD.tmp$Pseudotime,5/6))) %>% summarise(EMT_score_new=mean(EMT_score))
  EMT_compare_E=rbind(EMT_compare_E, pD.E)
  EMT_compare_H1=rbind(EMT_compare_H1, pD.H1)
  EMT_compare_H2=rbind(EMT_compare_H2, pD.H2)
  EMT_compare_H3=rbind(EMT_compare_H3, pD.H3)
  EMT_compare_H4=rbind(EMT_compare_H4, pD.H4)
  EMT_compare_M=rbind(EMT_compare_M, pD.M)
}


EMT_compare=data.frame(
  EMT_score_new=c(EMT_compare_E$EMT_score_new,EMT_compare_H1$EMT_score_new,EMT_compare_H2$EMT_score_new,EMT_compare_H3$EMT_score_new,EMT_compare_H4$EMT_score_new,EMT_compare_M$EMT_score_new),
  Group=c(rep('E',13),rep('H1',13),rep('H2',13),rep('H3',13),rep('H4',13),rep('M',13)))


EMT_compare$Phase[EMT_compare$Group=='E']=1
EMT_compare$Phase[EMT_compare$Group=='H1']=2
EMT_compare$Phase[EMT_compare$Group=='H2']=3
EMT_compare$Phase[EMT_compare$Group=='H3']=4
EMT_compare$Phase[EMT_compare$Group=='H4']=5
EMT_compare$Phase[EMT_compare$Group=='M']=6

ggplot(EMT_compare, aes(Phase, EMT_score_new, color=as.factor(Phase)))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  geom_jitter(position=position_jitter(width = 0.3, height=0),alpha=1,size=1.5)+
  geom_smooth(method = "lm", se=T, color="black", formula = y~x, size=0.75)+
  ggpubr::stat_cor(data=EMT_compare, mapping = aes(Phase,EMT_score_new), method='pearson', size=3, inherit.aes = F)+
  theme_bw()+
  scale_color_d3("category20")+
  scale_x_continuous(breaks = seq(1,6, by=1))+
  theme(axis.text=element_text(face = "bold"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+
  labs(title = "EMT score from E to M phase")
ggsave("palantir_result.pdf",width = 5,height = 3.8,units = "in")




###CellRank#######
Pt01_cr=read.csv("Pt01_cellrank.csv")
Pt02_cr=read.csv("Pt02_cellrank.csv")
Pt03_cr=read.csv("Pt03_cellrank.csv")
Pt04_cr=read.csv("Pt04_cellrank.csv")
Pt05_cr=read.csv("Pt05_cellrank.csv")
Pt06_cr=read.csv("Pt06_cellrank.csv")
Pt07_cr=read.csv("Pt07_cellrank.csv")
Pt08_cr=read.csv("Pt08_cellrank.csv")
Pt10_cr=read.csv("Pt10_cellrank.csv")
Pt11_cr=read.csv("Pt11_cellrank.csv")
Pt15_cr=read.csv("Pt15_cellrank.csv")
Pt17_cr=read.csv("Pt17_cellrank.csv")
Pt18_cr=read.csv("Pt18_cellrank.csv")


cr=Pt01_cr %>% rbind(Pt02_cr) %>% rbind(Pt03_cr) %>% rbind(Pt04_cr) %>% rbind(Pt05_cr) %>% 
  rbind(Pt06_cr) %>% rbind(Pt07_cr) %>% rbind(Pt08_cr) %>% rbind(Pt10_cr) %>% rbind(Pt11_cr) %>% 
  rbind(Pt15_cr) %>% rbind(Pt17_cr) %>% rbind(Pt18_cr)


head(cr)
cr_palantir=cr %>% group_by(patient,lesions) %>% summarise(cellrank=mean(term_states_fwd_probs, na.rm = TRUE))
cr_palantir<-cr_palantir %>% arrange(lesions,patient) 

wilcox.test(cr_palantir$cellrank[1:13],cr_palantir$cellrank[14:26],paired = T,alternative = c("less")) 

ggpaired(cr_palantir, x = 'lesions', y = 'cellrank',
         color='lesions',fill = 'lesions', line.color = "#7f7f7f",linetype ='dashed', line.size = 0.4, point.size = 2)+   
  scale_color_igv()+
  scale_fill_manual(values = c('#5050FF80','#CE3D3280'))+
  labs(x = '', y = 'CellRank pseudotime', title = '', subtitle = '')+geom_signif(xmin=1,
                                                                                 xmax=2, 
                                                                                 y_position=1.05,
                                                                                 annotations = 'p=0.029')
ggsave("CellRank_pseudotime_PT_vs_LM.pdf",width = 2.5,height = 4,units = "in")




pData.new=as.data.frame(matrix(nrow=0, ncol = 7))
colnames(pData.new)=c(colnames(cr),'Probability')
pt.list=unique(cr$patient)
for (i in pt.list) {
  pD.tmp=subset(cr, patient %in% i)
  pD.tmp$Probability=(pD.tmp$term_states_fwd_probs-min(pD.tmp$term_states_fwd_probs))/(max(pD.tmp$term_states_fwd_probs)-min(pD.tmp$term_states_fwd_probs))*100 
  pData.new=rbind(pData.new, pD.tmp) 
}

pData.new$phase[pData.new$Probability<=100/6]='E'
pData.new$phase[pData.new$Probability>100/6 & pData.new$Probability<=200/6]='H1'
pData.new$phase[pData.new$Probability>200/6 & pData.new$Probability<=300/6]='H2'
pData.new$phase[pData.new$Probability>300/6 & pData.new$Probability<=400/6]='H3'
pData.new$phase[pData.new$Probability>400/6 & pData.new$Probability<=500/6]='H4'
pData.new$phase[pData.new$Probability>500/6 & pData.new$Probability<=100]='M'


EMT_compare_E=c()
EMT_compare_H1=c()
EMT_compare_H2=c()
EMT_compare_H3=c()
EMT_compare_H4=c()
EMT_compare_M=c()

pt.list=unique(cr$patient)
for (i in pt.list) {
  pD.tmp=subset(pData.new, patient%in%i)
  pD.E =pD.tmp%>% subset(phase=='E') %>% summarise(EMT_score_new=mean(EMT_score[EMT_score<quantile(EMT_score,1/6)]))
  pD.H1 =pD.tmp%>% subset(phase=='H1')  %>% summarise(EMT_score_new=mean(EMT_score[EMT_score<quantile(EMT_score,2/6)]))
  pD.H2 =pD.tmp%>% subset(phase=='H2') %>% summarise(EMT_score_new=mean(EMT_score[EMT_score<quantile(EMT_score,3/6)]))
  pD.H3 =pD.tmp%>% subset(phase=='H3') %>% summarise(EMT_score_new=mean(EMT_score[EMT_score<quantile(EMT_score,4/6)]))
  pD.H4 =pD.tmp%>% subset(phase=='H4') %>% summarise(EMT_score_new=mean(EMT_score[EMT_score<quantile(EMT_score,5/6)]))
  pD.M =pD.tmp%>% subset(phase=='M') %>% summarise(EMT_score_new=mean(EMT_score[EMT_score<quantile(EMT_score,6/6)]))
  EMT_compare_E=rbind(EMT_compare_E, pD.E)
  EMT_compare_H1=rbind(EMT_compare_H1, pD.H1)
  EMT_compare_H2=rbind(EMT_compare_H2, pD.H2)
  EMT_compare_H3=rbind(EMT_compare_H3, pD.H3)
  EMT_compare_H4=rbind(EMT_compare_H4, pD.H4)
  EMT_compare_M=rbind(EMT_compare_M, pD.M)
}


EMT_compare=data.frame(
  EMT_score_new=c(EMT_compare_E$EMT_score_new,EMT_compare_H1$EMT_score_new,EMT_compare_H2$EMT_score_new,EMT_compare_H3$EMT_score_new,EMT_compare_H4$EMT_score_new,EMT_compare_M$EMT_score_new),
  Group=c(rep('E',13),rep('H1',13),rep('H2',13),rep('H3',13),rep('H4',13),rep('M',13)))


EMT_compare$Phase[EMT_compare$Group=='E']=1
EMT_compare$Phase[EMT_compare$Group=='H1']=2
EMT_compare$Phase[EMT_compare$Group=='H2']=3
EMT_compare$Phase[EMT_compare$Group=='H3']=4
EMT_compare$Phase[EMT_compare$Group=='H4']=5
EMT_compare$Phase[EMT_compare$Group=='M']=6


ggplot(EMT_compare, aes(Phase, EMT_score_new, color=as.factor(Phase)))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  geom_jitter(position=position_jitter(width = 0.3, height=0),alpha=1,size=1.5)+
  geom_smooth(method = "lm", se=T, color="black", formula = y~x, size=0.75)+
  ggpubr::stat_cor(data=EMT_compare, mapping = aes(Phase,EMT_score_new), method='pearson', size=3, inherit.aes = F)+
  theme_bw()+
  scale_color_d3("category20")+
  scale_x_continuous(breaks = seq(1,6, by=1))+
  theme(axis.text=element_text(face = "bold"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+
  labs(title = "EMT score from E to M phase")
ggsave("cellrank_result.pdf",width = 5,height = 3.8,units = "in")



