library(Seurat)
library(ggsci)
library(dplyr)
library(ggplot2)


Myeloid.major=readRDS("/Users/a1502/Desktop/Liuxin/001-Data/009-analysis-outs/CRC60samples/18_patients_filter/Myeloid_filter/myeloid_major_20220104_final.rds")
Myeloid.major=readRDS("/Users/xinliu/Desktop/Liuxin/001-Data/009-analysis-outs/CRC60samples/18_patients_filter/Myeloid_filter/myeloid_major_20220104_final.rds")



######Fig S9C barplot######
ratio_data=as.data.frame(table(Myeloid.major$cluster_name, Myeloid.major$lesions))
colnames(ratio_data)=c('cluster','lesions','Freq')
ratio_data$cluster=factor(ratio_data$cluster,levels=c("Macro_c1_CCL18","Macro_c3_CD163L1","Macro_c5_ANGPTL4","Macro_c7_CXCL10","Macro_c10_MT1G","Macro_c11_HSPH1","Macro_c12_CD5L",
                                                      "Mono_c2_FCN1","cDC2_c4_CLEC10A","DC_c8_LAMP3","pDC_c9_LILRA4","cDC1_c13_CLEC9A","Cycling_Myeloid_c6"),ordered = T)

ggplot(ratio_data, aes(fill=lesions, y=cluster, x=Freq)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_d3('category20')+theme_light()+
  theme(axis.title =element_text(size = 8),axis.text =element_text(size = 8, color = 'black'))+
  theme(axis.text.x = element_blank(), strip.background = element_blank(),strip.text = element_blank(),
        axis.ticks.x = element_blank())+scale_x_reverse()+
  xlab('Proportion')+ylab('')
ggsave("Fig S9C.pdf",width = 5,height = 4,units = "in")







#####Fig 8D######
total_meta = Myeloid.major@meta.data
head(total_meta)
unique(total_meta$lesions)
unique(total_meta$location)

cluster_name = levels(factor(total_meta$cluster_name))

orig.ident<-levels(factor(total_meta$orig.ident))
sample_list<-unique(orig.ident)

Myeloid_data_prop<-tibble(cluster=cluster_name)

for (i in 1:length(sample_list)) {
  tmp1<-as.data.frame(table(total_meta$cluster_name[total_meta$orig.ident==sample_list[i]]))
  colnames(tmp1)<-c("cluster",sample_list[i])
  Myeloid_data_prop<-Myeloid_data_prop%>%left_join(tmp1,by="cluster")
}
Myeloid_data_prop<-as.data.frame(Myeloid_data_prop)
Myeloid_data_prop[is.na(Myeloid_data_prop)]<-0
for (i in 2:(length(sample_list)+1)) {
  Myeloid_data_prop[,i]<-Myeloid_data_prop[,i]/sum(Myeloid_data_prop[,i])
}
colnames(Myeloid_data_prop)

Myeloid_data_prop<-Myeloid_data_prop%>%
  tidyr::pivot_longer(
    cols = "M5208-1":"X5226-5",
    names_to = "sample_ID",
    values_to = "cell_prop"
  )

Myeloid_data_prop$lesions<-"Colon_tumor"
Myeloid_data_prop$lesions[Myeloid_data_prop$sample_ID%in%unique(Myeloid.major@meta.data$orig.ident[Myeloid.major@meta.data$lesions=="Liver_meta"])]<-"Liver_meta"
Myeloid_data_prop$lesions<-factor(Myeloid_data_prop$lesions,level=c("Colon_tumor","Liver_meta"))
unique(Myeloid_data_prop$cluster)
Myeloid_data_prop=filter(Myeloid_data_prop, cluster%in%c('Macro_c1_CCL18','Macro_c3_CD163L1','Macro_c7_CXCL10',
                                                         'Macro_c10_MT1G','Macro_c11_HSPH1','Macro_c12_CD5L'))

Myeloid_data_prop$cluster=factor(Myeloid_data_prop$cluster,level=c(unique(Myeloid_data_prop$cluster)))


library(ggpubr)
Myeloid_data_prop = Myeloid_data_prop%>%mutate(a=as.numeric(cluster), 
                                           b=ifelse(lesions=='Colon_tumor', -0.2, 0.2)) 

ggplot(Myeloid_data_prop,aes(x=cluster,y=cell_prop,color=lesions))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(aes(a+b, cell_prop, color=lesions),position=position_jitter(width = 0.1, height=0),alpha=0.5,size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,vjust=0.9,hjust=0.9),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black")) + 
  stat_compare_means(method = "wilcox.test",size=2.5,label = "p.format",angle=30,vjust = 2)+   
  labs(title = "")+ylab('Cell proportion in Myeloid cells')+xlab('')

ggsave("Fig 8D.pdf",width = 4,height = 3,units = "in")


