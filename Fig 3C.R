library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsci)

Tcell=readRDS("Tcell_in_ann_new_1227.rds")

total_meta = Tcell@meta.data
cluster_name = levels(factor(total_meta$cluster_name))

a <- total_meta %>% group_by(lesions,cluster_name) %>%summarise(a=n())
ac <- total_meta %>% group_by(lesions) %>% summarise(ac=n())
ab <- total_meta %>% group_by(cluster_name) %>% summarise(ab=n())
abcd = dim(total_meta)[1]
dist.data <- merge(a, ac)
dist.data <- merge(dist.data, ab)

dist.data['b'] = dist.data$ab - dist.data$a       
dist.data['c'] = dist.data$ac - dist.data$a
dist.data['d'] = abcd - dist.data$ab - dist.data$c

head(dist.data)

# chisq.test
x=as.numeric(dist.data[1, c('a', 'c', 'b', 'd')])
x = chisq.test(matrix(x, nrow=2,ncol=2))

dist.data[,c('p.value','Ea', 'Ec', 'Eb', 'Ed')] = t(apply(dist.data, 1, function(x){
  x=as.numeric(x[c('a', 'c', 'b', 'd')])
  k.test = chisq.test(matrix(x, nrow=2,ncol=2))
  return(c(k.test$p.value, as.vector(k.test$expected)))
}))
dist.data['Reo'] = dist.data$a / dist.data$Ea
library(ggplot2)
library(reshape2)
library(pheatmap)
dist.heatmap.data <- dist.data[,c('cluster_name', 'lesions', 'Reo')]

dist.heatmap<-data.frame(matrix(data=0, nrow = length(cluster_name), ncol =2))
rownames(dist.heatmap) = cluster_name
colnames(dist.heatmap) = c('Colon_tumor','Liver_meta')
for(i in 1:dim(dist.heatmap.data)[1]){
  ct = as.character(dist.heatmap.data[i, 'cluster_name'])
  ts = dist.heatmap.data[i, 'lesions']
  reo = dist.heatmap.data[i, 'Reo']
  dist.heatmap[ct, ts] = reo
}
dist.heatmap1<-dist.heatmap
dist.heatmap1=dist.heatmap1[sapply(cluster_name,function(e){which(rownames(dist.heatmap1)==e)}),]


annote.heatmap = dist.heatmap
annote.heatmap[annote.heatmap>1] = '+++'
annote.heatmap[annote.heatmap<=1&annote.heatmap>0.8] = '++'
annote.heatmap[annote.heatmap<=0.8&annote.heatmap>=0.2] = '+'
annote.heatmap[annote.heatmap<0.2&annote.heatmap>0] = '+/-'
annote.heatmap[annote.heatmap==0] = '-'


dist.heatmap1=dist.heatmap1 %>% arrange(desc(Colon_tumor)) 
cell_cluster_dist=pheatmap(dist.heatmap1,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           treeheight_row = 0,
                           color = colorRampPalette(c("#FFFFFC","#FF9900"))(50),
                           display_numbers =round(as.matrix(dist.heatmap1),2),
                           cellwidth = 30, cellheight = 16, fontsize = 12,
                           border_color = '#ffffff',
                           angle_col='45',
                           main = '')

cell_cluster_dist


#####convert
a <- total_meta %>% group_by(lesions,cluster_name) %>%summarise(a=n()) 
ac <- total_meta %>% group_by(lesions) %>% summarise(ac=n())  
ab <- total_meta %>% group_by(cluster_name) %>% summarise(ab=n()) 
abcd = dim(total_meta)[1]  
dist.data <- merge(a, ac)
dist.data <- merge(dist.data, ab)

dist.data['b'] = dist.data$ab - dist.data$a       
dist.data['c'] = dist.data$ac - dist.data$a
dist.data['d'] = abcd - dist.data$ab - dist.data$c

head(dist.data)

dist.data[,c('p.value','Ea', 'Ec', 'Eb', 'Ed')] = t(apply(dist.data, 1, function(x){
  x=as.numeric(x[c('a', 'c', 'b', 'd')])
  k.test = chisq.test(matrix(x, nrow=2,ncol=2))
  return(c(k.test$p.value, as.vector(k.test$expected)))
}))
dist.data['Reo'] = dist.data$a / dist.data$Ea

library(ggplot2)
library(reshape2)
library(pheatmap)
head(dist.data)
dist.heatmap.data <- dist.data[,c('cluster_name', 'lesions','p.value','Reo')]
head(dist.heatmap.data)

dist.heatmap<-data.frame(matrix(data=0, nrow = length(cluster_name), ncol =3))
rownames(dist.heatmap) = cluster_name
colnames(dist.heatmap) = c('Colon_tumor','Liver_meta','p.value')
rownames(dist.heatmap.data)<-NULL
for(i in 1:dim(dist.heatmap.data)[1]){
  ct = as.character(dist.heatmap.data[i, 'cluster_name'])
  ts = dist.heatmap.data[i, 'lesions']
  reo = dist.heatmap.data[i, 'Reo']
  p=dist.heatmap.data[i, 'p.value']
  dist.heatmap[ct, as.vector(ts)] = reo
  dist.heatmap[ct, 'p.value'] = p
}

dist.heatmap1<-dist.heatmap[,1:2]
dist.heatmap1<-dist.heatmap1 %>% arrange(desc(Colon_tumor)) 


annote.heatmap = dist.heatmap
annote.heatmap$annote[annote.heatmap$p.value>=0.01] = ''
annote.heatmap$annote[annote.heatmap$p.value<0.05] = '*'
annote.heatmap$annote[annote.heatmap$p.value<0.01] = '**'
annote.heatmap$annote[annote.heatmap$p.value<0.001] = '***'
annote.heatmap<-annote.heatmap %>% arrange(desc(Colon_tumor)) 

annote.heatmap$Colon_tumor[annote.heatmap$Colon_tumor<1] = ""
annote.heatmap$Colon_tumor[annote.heatmap$Colon_tumor>1] = annote.heatmap$annote[annote.heatmap$Colon_tumor>1]
annote.heatmap$Liver_meta[annote.heatmap$Liver_meta<1] = ""
annote.heatmap$Liver_meta[annote.heatmap$Liver_meta>1] = annote.heatmap$annote[annote.heatmap$Liver_meta>1]
annote.heatmap$Colon_tumor[annote.heatmap$Colon_tumor<annote.heatmap$Liver_meta] = ""
annote.heatmap$Liver_meta[annote.heatmap$Liver_meta<annote.heatmap$Colon_tumor] = ""

annote.heatmap<-annote.heatmap[,c('Colon_tumor','Liver_meta')]


cell_cluster_dist=pheatmap(dist.heatmap1,
                           cluster_rows = FALSE,
                           cluster_cols = FALSE,
                           treeheight_row = 0,
                           color = colorRampPalette(c("#FFFFFC","#FF9900"))(50),
                           display_numbers =as.matrix(annote.heatmap),
                           cellwidth = 30, cellheight = 16, fontsize = 12,
                           border_color = '#ffffff',
                           angle_col='45',
                           main = '')
cell_cluster_dist



#####fig 3C######
total_meta = Tcell@meta.data
head(total_meta)
unique(total_meta$lesions)
unique(total_meta$location)

cluster_name = levels(factor(total_meta$cluster_name))

orig.ident<-levels(factor(total_meta$orig.ident))
sample_list<-unique(orig.ident)

Tcell_data_prop<-tibble(cluster=cluster_name)

for (i in 1:length(sample_list)) {
  tmp1<-as.data.frame(table(total_meta$cluster_name[total_meta$orig.ident==sample_list[i]]))
  colnames(tmp1)<-c("cluster",sample_list[i])
  Tcell_data_prop<-Tcell_data_prop%>%left_join(tmp1,by="cluster")
}
Tcell_data_prop<-as.data.frame(Tcell_data_prop)
Tcell_data_prop[is.na(Tcell_data_prop)]<-0
for (i in 2:(length(sample_list)+1)) {
  Tcell_data_prop[,i]<-Tcell_data_prop[,i]/sum(Tcell_data_prop[,i])
}
colnames(Tcell_data_prop)

Tcell_data_prop<-Tcell_data_prop%>%
  tidyr::pivot_longer(
    cols = "M5208-1":"X5226-5",
    names_to = "sample_ID",
    values_to = "cell_prop"
  )

Tcell_data_prop$lesions<-"Colon_tumor"
Tcell_data_prop$lesions[Tcell_data_prop$sample_ID%in%unique(Tcell@meta.data$orig.ident[Tcell@meta.data$lesions=="Liver_meta"])]<-"Liver_meta"
Tcell_data_prop$lesions<-factor(Tcell_data_prop$lesions,level=c("Colon_tumor","Liver_meta"))
Tcell_data_prop$cluster=factor(Tcell_data_prop$cluster,level=c('CD4_IL17A_Th17','CD4_CCR7_Tn-like1','CD8_SLC4A10_MAIT',
                                                               'CD8_CXCL13_Tex','TYROBP_NK&NKT','CD8_GZMK_Tem',
                                                               'Cycling_T','CD4_IL7R_Tn-like2','CD4_FOXP3_Treg'))

library(ggpubr)
Tcell_data_prop = Tcell_data_prop%>%mutate(a=as.numeric(cluster), 
                                           b=ifelse(lesions=='Colon_tumor', -0.2, 0.2))  

ggplot(Tcell_data_prop,aes(x=cluster,y=cell_prop,color=lesions))+
  #stat_boxplot(geom = "errorbar")+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_d3("category20")+
  geom_jitter(aes(a+b, cell_prop, color=lesions),position=position_jitter(width = 0.1, height=0),alpha=0.5,size=0.8)+
  theme_bw()+
  theme(axis.text=element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,vjust=0.9,hjust=0.9),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black")) + 
  stat_compare_means(method = "wilcox.test",size=2.5,label = "p.format",angle=30,vjust = 2)+   
  labs(title = "")+ylab('Cell proportion in T cells')+xlab('')

ggsave("fig 3C.pdf",width = 6,height = 5,units = "in")


#RORC UMAP for revise######
CD4T=Tcell %>% subset(cycling_sub%in%c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                       'CD4_Cycling_T'))

FeaturePlot(CD4T, features=c('RORC'))

CD4T$Th17='Not Th17'
CD4T$Th17[CD4T$cluster_name=='CD4_IL17A_Th17']='Th17'

a=FeaturePlot(CD4T, features=c('RORC'))
b=DimPlot(CD4T,group.by = 'Th17')+scale_color_manual(values = c('gray','red'))+NoLegend()

cowplot::plot_grid(a, b, ncol=2)

ggsave("CD4T_RORC_umap.pdf",width = 6,height = 3,units = "in")



#####RORC expression for revise######

Th17=CD4T %>% subset(cycling_sub%in%c('CD4_IL17A_Th17'))
CD4Th17_RORC=Th17 %>% subset(RORC>0)
dim(CD4Th17_RORC)[2]/dim(Th17)[2]
#0.2211862 of Th17 express RORC


library(tidyr)
library(cowplot)
avg_expression <- Tcell@assays$RNA@data[c('RORC','RORA'), ] %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = Tcell$cycling_sub) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean))

avg_expression_long <- avg_expression %>%
  gather(gene, mean, -cluster)

exp_data <- Tcell@assays$RNA@data[c('RORC','RORA'), ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = Tcell$cycling_sub) %>%
  gather(gene, exp, -cluster)

final_data <- left_join(avg_expression_long, exp_data, by = c("cluster", "gene"))

final_data=final_data %>% filter(cluster%in%c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                              'CD4_Cycling_T'))
final_data$cluster=factor(final_data$cluster,levels = c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                                        'CD4_Cycling_T'))

ggplot(final_data, aes(cluster, exp,fill = mean)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE,color="white")+
  scale_fill_gradient(low = "#f2e80d",high = "#e65224") +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(gene), scales = "free", switch = "y") +
  theme_cowplot(font_size = 10) +
  theme(panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "plain"),
        strip.text.y.left = element_text(angle = 0,hjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(size=0.1),
        axis.ticks = element_line(size=0.1)) +
  labs(title = '')+ theme(plot.title = element_text(hjust = 0.5))+xlab('')

ggsave("RORC_RORA_CD4_Tcell.pdf",width = 3,height = 2.5,units = "in")


#####CD4 IL17 gene family for revise######
CD4_T=Tcell %>% subset(cycling_sub%in%c('CD4_CCR7_Tn-like1',"CD4_IL7R_Tn-like2",'CD4_FOXP3_Treg','CD4_IL17A_Th17'))  
genelist=c('IL17A','IL17B','IL17C','IL17D','IL17F') 
library(tidyr)
library(cowplot)
avg_expression <- CD4_T@assays$RNA@data[genelist, ] %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = CD4_T$cycling_sub) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean))

avg_expression_long <- avg_expression %>%
  gather(gene, mean, -cluster)

exp_data <- CD4_T@assays$RNA@data[genelist, ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = CD4_T$cycling_sub) %>%
  gather(gene, exp, -cluster)

final_data <- left_join(avg_expression_long, exp_data, by = c("cluster", "gene"))

final_data$cluster=factor(final_data$cluster,levels = c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2'))


ggplot(final_data, aes(cluster, exp,fill = mean)) +
  geom_violin(scale = "width", adjust = 2, trim = TRUE,color="white")+
  scale_fill_gradient(low = "#f2e80d",high = "#e65224") +
  scale_y_continuous(expand = c(0, 0), position="right")+ #, labels = function(x)c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(gene), scales = "free", switch = "y") +
  theme_cowplot(font_size = 8) +
  theme(panel.spacing = unit(0, "lines"),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "plain"),
        strip.text.y.left = element_text(angle = 0,hjust = 1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(size=0.1),
        axis.ticks = element_line(size=0.1)) +
  labs(title = '')+ theme(plot.title = element_text(hjust = 0.5))+xlab('')

ggsave("CD4_Tcell_IL17ABCDF_expression.pdf",width = 3,height = 3,units = "in")



#Tcell_Mito_percent for revise#####
VlnPlot(Tcell, features=c('percent.mt'), group.by = 'cluster_name',pt.size=0)
ggsave("Tcell_Mito_percent.pdf",width = 6,height = 4,units = "in")


quantile(Tcell$percent.mt, probs = 0.9)
#10.4093 
#about 90% of T cells had mitochondrial gene content below 10%. 






