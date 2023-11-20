library(ggplot2)
library(dplyr)
library(ggsci)
library(ggpubr)
library(patchwork)
df=read.csv('FAH_SYSU_NPX.csv')
sample_use=unique(df$SampleID) %>% setdiff(c('R21069880','R21069884','R21069885','R21069886','R21069887','R21069889'))
df=subset(df, df$SampleID %in% sample_use)

df$Group='HCC'
df$Group[df$SampleID%in%c("test_1","test_9")]='Normal'
df$Group[df$SampleID%in%c('R21069871','R21069872','R21069873','R21069874','R21069875','R21069876','R21069877','R21069878','R21069879','R21069880',
                          'R21069881','R21069882','R21069883','R21069884','R21069885','R21069886','R21069887','R21069889')]='CRLM'
df$Group=factor(df$Group,levels=c('Normal','HCC','CRLM'),ordered = T)

compaired=list(c('CRLM','Normal'))

#Fig 8L####
ggplot(df %>% filter(Assay=='SPP1')%>%filter(Group%in%c('Normal','CRLM')),aes(x=Group,y=NPX,color=Group))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  scale_color_igv()+
  geom_jitter(position=position_jitter(width = 0.3, height=0),alpha=1,size=2)+
  theme_bw()+
  theme(axis.text=element_text(colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0, hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,colour = "black"))+
  labs(title = "SPP1")+xlab('')+ylab('NPX')+ylim(c(0,3.25))+
  geom_signif(comparisons = compaired, tip_length = 0.01,step_increase = 0.01, map_signif_level = F, test = 'wilcox.test',color='black')
ggsave("Fig 8L.pdf",width = 4,height = 6,units = "in")


