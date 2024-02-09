library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)
library(psych)
library(qgraph)
library(igraph)
library(tidyverse)

mypvals_PT <- read.delim(paste0("PT_pvalues.txt"), check.names = FALSE)
mymeans_PT <- read.delim(paste0("PT_means.txt"), check.names = FALSE)
mysigmeans_PT <- read.delim(paste0("PT_significant_means.txt"), check.names = FALSE)

mypvals_LM <- read.delim(paste0("LM_pvalues.txt"), check.names = FALSE)
mymeans_LM <- read.delim(paste0("LM_means.txt"), check.names = FALSE)
mysigmeans_LM <- read.delim(paste0("LM_significant_means.txt"), check.names = FALSE)

myfocus=c('CCL4_SLC7A1','CCL3L1_DPP4','CXCL2_DPP4','SPP1_CD44','CCL3_IDE','CXCL10_DPP4','GRN_SORT1') 

mysigmeans_PT %>% dplyr::filter(interacting_pair %in% c(myfocus))%>%
  dplyr::select("interacting_pair","Macro_c3_CD163L1|Tumor","Macro_c1_CCL18|Tumor","Macro_c7_CXCL10|Tumor") %>%  
  reshape2::melt() -> meansdf
colnames(meansdf)<- c("interacting_pair","CC","means")

mypvals_PT %>% dplyr::filter(interacting_pair %in% c(myfocus))%>%
  dplyr::select("interacting_pair","Macro_c3_CD163L1|Tumor","Macro_c1_CCL18|Tumor","Macro_c7_CXCL10|Tumor")%>%  
  reshape2::melt()-> pvalsdf
colnames(pvalsdf)<- c("interacting_pair","CC","pvals")

pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf_PT <- merge(pvalsdf,meansdf,by = "joinlab")
pldf_PT$CC.x=paste0(pldf_PT$CC.x,' (PT)')



mysigmeans_LM %>% dplyr::filter(interacting_pair %in% c(myfocus))%>%
  dplyr::select("interacting_pair","Macro_c3_CD163L1|Tumor","Macro_c1_CCL18|Tumor","Macro_c7_CXCL10|Tumor") %>%  
  reshape2::melt() -> meansdf
colnames(meansdf)<- c("interacting_pair","CC","means")

mypvals_LM %>% dplyr::filter(interacting_pair %in% c(myfocus))%>%
  dplyr::select("interacting_pair","Macro_c3_CD163L1|Tumor","Macro_c1_CCL18|Tumor","Macro_c7_CXCL10|Tumor")%>%  
  reshape2::melt()-> pvalsdf
colnames(pvalsdf)<- c("interacting_pair","CC","pvals")

pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf_LM <- merge(pvalsdf,meansdf,by = "joinlab")
pldf_LM$CC.x=paste0(pldf_LM$CC.x,' (LM)')



pldf_total=rbind(pldf_PT,pldf_LM)

pldf_total$CC.x=factor(pldf_total$CC.x,levels = c('Macro_c1_CCL18|Tumor (PT)',
                                                  'Macro_c7_CXCL10|Tumor (PT)',
                                                  'Macro_c3_CD163L1|Tumor (PT)',
                                                  'Macro_c1_CCL18|Tumor (LM)',
                                                  'Macro_c7_CXCL10|Tumor (LM)',
                                                  'Macro_c3_CD163L1|Tumor (LM)'),ordered = T)
pldf_total[is.na(pldf_total)]=0
pldf_total$pvals[pldf_total$pvals>0.05]=1
pldf_total$Interaction_score=pldf_total$means*4/max(pldf_total$means)


##Fig S9D######
pldf_total %>% 
  ggplot(aes(CC.x,interacting_pair.x))+ 
  geom_point(aes(color=Interaction_score,size=-log10(pvals+0.0001))) +
  #scale_size_continuous(range = c(1,6))+
  scale_size_continuous(breaks=c(0, -log10(0.05+0.0001), -log10(0.01+0.0001), -log10(0.001+0.0001)), labels=c('ns','*', '**', '***')) +
  guides(size=guide_legend(title='Significance'))+
  scale_color_gradient2(high="red",low ="white",mid='yellow',midpoint = 2)+ theme_bw()+ theme(panel.grid.major = element_blank())+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.9,vjust = 0.9,face = 'bold',color = 'black'),
        axis.text.y = element_text(face = 'bold',color = 'black')) +
  geom_vline(xintercept = c(3.5),size=0.1,linetype=2)+
  ylab('interacting_pair')+xlab('')

ggsave("Fig S9D.pdf",width = 5,height = 6,units = "in")
















