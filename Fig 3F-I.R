library(dplyr)
library(patchwork)
library(Seurat)
library(ggplot2)
library(future) # multiprocess
library(cowplot)
library(ggsci)
library(gg.gap)
library(igraph)

#####Fig 3G#######
Tcell_srt = readRDS("Tcell_in_ann_new_1227_withTCR.rds")
library(ComplexHeatmap)
library(reshape2)

metadata = subset(Tcell_srt@meta.data, enrich_singlet=='enriched')
metadata = subset(metadata, cluster_name%in%c( 'CD4_IL17A_Th17','CD4_FOXP3_Treg'))
metadata$group = factor(metadata$cluster_name, levels = c( 'CD4_IL17A_Th17','CD4_FOXP3_Treg'))

share_mat = metadata%>%
  group_by(group, raw_clonotype_id) %>%
  summarise(cell_num=n()) %>% as.data.frame()
share_mat = as.data.frame(dcast(share_mat, formula = raw_clonotype_id~group))
rownames(share_mat) = as.character(share_mat$raw_clonotype_id)
share_mat = share_mat[,-1]
share_mat[is.na(share_mat)] = 0
share_mat[share_mat>5] = 5
share_mat = share_mat[rowSums(share_mat)!=0,]
colnames(share_mat) = c('left', 'right')
##
left_unique = share_mat[share_mat$left>0 & share_mat$right==0, ] %>% arrange(-left, right)
right_unique = share_mat[share_mat$left==0 & share_mat$right>0, ] %>% arrange(-left, right)
share_unique = share_mat[share_mat$left>0 & share_mat$right>0, ] %>% arrange(-left, -right)

new_share_mat = rbind(share_unique, left_unique, right_unique)

cn_col = structure(c('gray',
                     '#FFEB97', 
                     '#FCCC00', 
                     '#ec9336',
                     '#FC0000',
                     '#7d170e'), names=0:5)
Heatmap(new_share_mat,
        cluster_rows =F, cluster_columns = F,
        show_row_names = F, show_column_names = T,
        col = cn_col,
        column_labels = c('CD4_IL17A_Th17','CD4_FOXP3_Treg'),
        column_names_side = "top",
        column_names_rot = 0,
        column_names_gp = gpar(fontsize = 16), 
        column_names_centered = TRUE,
        #rect_gp = gpar(col='white'),
        heatmap_legend_param = list(
          at = 0:5,
          labels = c('0', '1', '2', '3', '4','>=5'),
          title = "Cell_num",
          legend_height = unit(4, "cm"),
          # title_position = "lefttop-rot",
          color_bar = "discrete"
        ))



### Fig 3F#########
metadata = Tcell_srt@meta.data

tmp_tab = table(na.omit(metadata$raw_clonotype_id))
tmp_names1 = names(tmp_tab[tmp_tab>1])
tmp_names2 = names(tmp_tab[tmp_tab==1])
enriched_bd =  rownames(subset(metadata, raw_clonotype_id%in%tmp_names1))
singlet_bd =  rownames(subset(metadata, raw_clonotype_id%in%tmp_names2))

Tcell_srt[['enriched']] = 0
Tcell_srt@meta.data[enriched_bd, 'enriched'] = 1
Tcell_srt[['singlet']] = 0
Tcell_srt@meta.data[singlet_bd, 'singlet'] = 1

enriched_data = Tcell_srt@meta.data
Tcell_srt[['enrich_singlet']] = 'unmap'
Tcell_srt@meta.data[Tcell_srt$singlet==1,'enrich_singlet'] = 'singlet'
Tcell_srt@meta.data[Tcell_srt$enriched==1,'enrich_singlet'] = 'enriched'
Tcell_srt[['barcode']] = colnames(Tcell_srt)
clono_cell_num = Tcell_srt@meta.data%>% 
  group_by(raw_clonotype_id)%>%
  dplyr::summarise(n=n()) %>%
  right_join(Tcell_srt@meta.data[,c('barcode', 'raw_clonotype_id')])%>%as.data.frame()
rownames(clono_cell_num) = clono_cell_num$barcode
clono_cell_num[is.na(clono_cell_num$raw_clonotype_id), 'n']=0
Tcell_srt[['clono_cell_num']] = cut(clono_cell_num[colnames(Tcell_srt), 'n'],
                                    breaks = c(0,1,2,5,8,Inf),right=F )
Tcell_df=Tcell_srt %>% subset(cluster_name%in%c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_Cycling_T'))
DimPlot(Tcell_df, group.by = 'clono_cell_num', cols = c('#dddddd',
                                                        '#FFEB97',
                                                        '#FCCC00', 
                                                        '#ec9336',
                                                        '#FC0000'))+
  labs(title = '')
ggsave("Fig 3F.pdf",width = 7,height = 6,units = "in")



####Fig 3H####
# install.packages("devtools")
# devtools::install_github("Japrin/STARTRAC")
library("Startrac")
metadata = Tcell_srt@meta.data[grep('CD4',Tcell_srt$cluster_name),]
indat = metadata[,c('barcode', 'raw_clonotype_id','enrich_singlet',
                    'patient', 'orig.ident', 'cluster_name', 'cluster_name','lesions')]
colnames(indat) = c("Cell_Name","clone.id","clone.status","patient","sampleType","stype","majorCluster","loc")

indat = subset(indat, clone.status!='unmap')
indat$stype = sapply(indat$stype, function(x)strsplit(x,'_')[[1]][1])
indat[indat$clone.status=='singlet', 'clone.status'] = 'NoClonal'
indat[indat$clone.status=='enriched', 'clone.status'] = 'Clonal'
out <- Startrac.run(indat, proj="CRC", cores=NULL,verbose=F)

library("ggpubr")
library("ggplot2")
library("tictoc")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library(plyr)
library(doParallel)
library(data.table)
dd = as.data.table(out@cluster.sig.data)[aid!=out@proj,][order(majorCluster),]
expa_data = dd[index=='expa',]
expa_data = expa_data[!is.na(expa_data$value), ]
cluster_order = expa_data %>% group_by(majorCluster) %>% 
  dplyr::summarise(res=median(value)) %>% arrange(-res)
expa_data$majorCluster = factor(expa_data$majorCluster, levels=cluster_order$majorCluster)
ggboxplot(expa_data,
          x="majorCluster",y="value",palette = "npg",
          color = "darkgreen", add = "point", outlier.colour=NULL) +
  stat_compare_means(comparisons = list(c('CD4_IL17A_Th17', 'CD4_FOXP3_Treg'),
                                        c('CD4_IL17A_Th17', 'CD4_IL7R_Tn-like2'),
                                        c('CD4_FOXP3_Treg', 'CD4_IL7R_Tn-like2')),method = 'wilcox.test')+
  labs(title = 'Expansion')+xlab('')+ylab('STARTRAC-expa')+
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'None',
        plot.title = element_text(hjust = 0.5)) 
ggsave("Fig 3H.pdf",width = 4,height = 6,units = "in")



#####Fig 3I#####
dd = as.data.table(out@pIndex.sig.migr)[aid!=out@proj,][order(majorCluster),]
dd = dd[!is.na(dd$value), ]
cluster_order = dd %>% group_by(majorCluster) %>% 
  dplyr::summarise(res=median(value)) %>% arrange(-res)
dd$majorCluster = factor(dd$majorCluster, levels=cluster_order$majorCluster)

ggboxplot(dd,
          x="majorCluster",y="value",palette = "npg",
          color = "#56aaff", add = "point", outlier.colour=NULL) +
  stat_compare_means(comparisons = list(c('CD4_IL17A_Th17', 'CD4_FOXP3_Treg'),
                                        c('CD4_FOXP3_Treg', 'CD4_CCR7_Tn-like1'),
                                        c('CD4_IL17A_Th17', 'CD4_CCR7_Tn-like1')),
                     method = 'wilcox.test')+
  labs(title = 'Migration')+xlab('')+ylab('STARTRAC-migr')+
  theme(axis.text.x=element_text(angle = 60,hjust = 1),
        legend.position = 'None',
        plot.title = element_text(hjust = 0.5))     
ggsave("Fig 3I.pdf",width = 4,height = 6,units = "in")









