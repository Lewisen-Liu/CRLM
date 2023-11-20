library(Seurat)
library(dplyr)
library(ggsci)

Total<-readRDS("Total_20220309.rds")
table(Total$orig.ident,Total$patient)

######Fig S1B kBET#######
library(kBET)
table(Total$patient)
table(Total$Cluster)
Tcell<-Total %>% subset(Cluster%in%c('Tcell'))

data.PT<-GetAssayData(Tcell,slot='scale.data')
batch<-Tcell$patient

library(ggplot2)
library(ggpubr)
batch.estimate <- kBET(t(data.PT), batch, plot=F)
plot.data.PT1 <- data.frame(class=rep(c('observed', 'expected'), each=length(batch.estimate$stats$kBET.observed)), 
                            data.PT =  c(batch.estimate$stats$kBET.observed, batch.estimate$stats$kBET.expected))

Tcell_plot <- ggplot(plot.data.PT1, aes(class, data.PT,color=class)) + geom_boxplot() + scale_color_d3()+
  labs(x='Test', y='kBET rejection rate',title='T cell') +
  theme_bw() +  theme(axis.text = element_text(color='black'))+
  scale_y_continuous(limits=c(0,1.05))+ 
  stat_compare_means(method = "wilcox.test",size=4,label = "p.format", angle=0, vjust = 0.5, hjust=0)+NoLegend()
Tcell_plot
ggsave("Fig S1B Tcell_kbet.pdf",width = 2.5,height = 3,units = "in")



#Myeloid
Myeloid<-Total %>% subset(Cluster%in%c('Myeloid'))

data.Myeloid<-GetAssayData(Myeloid,slot='scale.data')
batch<-Myeloid$patient

library(ggplot2)
library(ggpubr)
batch.estimate.Myeloid <- kBET(t(data.Myeloid), batch, plot=F)
plot.data.Myeloid <- data.frame(class=rep(c('observed', 'expected'), each=length(batch.estimate.Myeloid$stats$kBET.observed)), 
                                data.Myeloid =  c(batch.estimate.Myeloid$stats$kBET.observed, batch.estimate.Myeloid$stats$kBET.expected))

Myeloid_plot <- ggplot(plot.data.Myeloid, aes(class, data.Myeloid ,color=class)) + geom_boxplot() + scale_color_d3()+
  labs(x='Test', y='kBET rejection rate',title='Myeloid') +
  theme_bw() +  theme(axis.text = element_text(color='black'))+
  scale_y_continuous(limits=c(0,1.05))+ 
  stat_compare_means(method = "wilcox.test",size=4,label = "p.format", angle=0, vjust = 0.5, hjust=0)+NoLegend()
Myeloid_plot
ggsave("Fig S1B Myeloid_kbet.pdf",width = 2.5,height = 3,units = "in")


#Bcell
Bcell<-Total %>% subset(Cluster%in%c('Bcell'))

data.Bcell<-GetAssayData(Bcell,slot='scale.data')
batch<-Bcell$patient

library(ggplot2)
library(ggpubr)
batch.estimate.Bcell <- kBET(t(data.Bcell), batch, plot=F)
plot.data.Bcell <- data.frame(class=rep(c('observed', 'expected'), each=length(batch.estimate.Bcell$stats$kBET.observed)), 
                                data.Bcell =  c(batch.estimate.Bcell$stats$kBET.observed, batch.estimate.Bcell$stats$kBET.expected))

Bcell_plot <- ggplot(plot.data.Bcell, aes(class, data.Bcell ,color=class)) + geom_boxplot() + scale_color_d3()+
  labs(x='Test', y='kBET rejection rate',title='B cell') +
  theme_bw() +  theme(axis.text = element_text(color='black'))+
  scale_y_continuous(limits=c(0,1.05))+ 
  stat_compare_means(method = "wilcox.test",size=4,label = "p.format", angle=0, vjust = 0.5, hjust=0)+NoLegend()
Bcell_plot
ggsave("Fig S1B Bcell_kbet.pdf",width = 2.5,height = 3,units = "in")



#Stromal
Stromal<-Total %>% subset(Cluster%in%c('Stromal'))

data.Stromal<-GetAssayData(Stromal,slot='scale.data')
batch<-Stromal$patient

library(ggplot2)
library(ggpubr)
batch.estimate.Stromal <- kBET(t(data.Stromal), batch, plot=F)
plot.data.Stromal <- data.frame(class=rep(c('observed', 'expected'), each=length(batch.estimate.Stromal$stats$kBET.observed)), 
                              data.Stromal =  c(batch.estimate.Stromal$stats$kBET.observed, batch.estimate.Stromal$stats$kBET.expected))

Stromal_plot <- ggplot(plot.data.Stromal, aes(class, data.Stromal ,color=class)) + geom_boxplot() + scale_color_d3()+
  labs(x='Test', y='kBET rejection rate',title='Stromal') +
  theme_bw() +  theme(axis.text = element_text(color='black'))+
  scale_y_continuous(limits=c(0,1.05))+ 
  stat_compare_means(method = "wilcox.test",size=4,label = "p.format", angle=0, vjust = 0.5, hjust=0)+NoLegend()
Stromal_plot
ggsave("Fig S1B Stromal_kbet.pdf",width = 2.5,height = 3,units = "in")




