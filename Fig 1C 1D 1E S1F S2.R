library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(dplyr)
library(ggsci)

process_seurat <-function(obj, dim=30, n.components=2,resolution=0.5, npcs=50){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = VariableFeatures(object = obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = npcs, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:dim)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:dim,reduction = "pca", n.components=n.components)
  return(obj)
}
regec_srt = readRDS("Regev_CRC_intratumor_immune.rds")
my_srt = readRDS("Total_20220309.rds")

my_srt = subset(my_srt, cluster_name!='Tumor cell')
my_srt = subset(my_srt, Cluster%in%c('Bcell', 'Tcell', 'Myeloid'))
regec_srt = process_seurat(regec_srt)

###fig 1C#####
Cluster_color=c('Tcell'='#b56cc1','Bcell'='#49a87e',
                'Myeloid'='#fec836', 'Stromal'='#ef5f5f')
a=DimPlot(my_srt, group.by='Cluster')+scale_color_manual(values = Cluster_color)
a
ggsave('fig 1C.pdf',
       a, width=220, height=220, units='mm', dpi = 600, bg = 'transparent')


### project UMAP
DefaultAssay(my_srt) = 'RNA'
regec_srt <- NormalizeData(regec_srt, verbose = FALSE)
regec_srt <- FindVariableFeatures(regec_srt, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

my_srt <- FindVariableFeatures(my_srt, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
anchors <- FindTransferAnchors(reference = my_srt, query = regec_srt, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = my_srt$cluster_name, dims = 1:30)
regec_srt <- AddMetaData(regec_srt, metadata = predictions)
# map umap
my_srt <- RunUMAP(my_srt, dims = 1:30, reduction = "pca", return.model = TRUE)
new_srt <- MapQuery(anchorset = anchors, 
                    reference = my_srt, 
                    query = regec_srt,
                    refdata = list(cluster_name = "cluster_name"), 
                    reference.reduction = "pca",
                    reduction.model = "umap",
                    projectumap.args=list(reduction.name='umap')
)


new_srt$sub_cluster = new_srt$predicted.cluster_name

project_srt <- merge(my_srt, new_srt, merge.dr = "umap")
Idents(project_srt) = project_srt$cluster_name
project_srt[['From']] = 'CRC'
project_srt@meta.data[colnames(regec_srt), 'From'] = 'Regev'

# saveRDS(mapping_srt.rds')
# project_srt=readRDS('mapping_srt.rds')

cell_mata = project_srt@meta.data %>%
  group_by(Cluster, cluster_name) %>%
  summarise(num=n()) %>% as.data.frame()
cluster_2_name = cell_mata$Cluster
names(cluster_2_name) = cell_mata$cluster_name

project_srt[['Group']] = project_srt$Cluster
project_srt@meta.data[is.na(project_srt$Group), 'Group'] = cluster_2_name[project_srt@meta.data[is.na(project_srt$Group),
                                                                                                'predicted.cluster_name']]
Idents(project_srt) = project_srt$Group
library(ggforce)

cell_color_liuxin=read.csv("color code.csv",header = F)

cell_color = c(cell_color_liuxin$V2, '#808080')
names(cell_color) = c(cell_color_liuxin$V1, 'CRC_cells')

project_srt[['CellType']] = project_srt$predicted.cluster_name
project_srt@meta.data[is.na(project_srt$predicted.cluster_name), 'CellType'] = 'CRC_cells'
cell_order = names(table(project_srt$predicted.cluster_name))
cell_order = c(cell_order, 'CRC_cells')

#fig 1D#####
a = DimPlot(project_srt, seed = 1, reduction = "umap", 
            group.by = "CellType",label = F, order = levels(my_srt$cluster_name), raster = F)+
  scale_color_manual(values = cell_color)+
  labs(title = 'Regev project')
a

ggsave('fig 1D.pdf',
       a, width=280, height=220, units='mm', dpi = 600, bg = 'transparent')



#fig 1E#####
cell_mata = project_srt@meta.data %>%
  group_by(Cluster, cluster_name) %>%
  summarise(num=n()) %>% as.data.frame()
cluster_2_name = cell_mata$Cluster
names(cluster_2_name) = cell_mata$cluster_name

project_srt[['Group']] = project_srt$Cluster
project_srt@meta.data[is.na(project_srt$Group), 'Group'] = cluster_2_name[project_srt@meta.data[is.na(project_srt$Group), 'predicted.cluster_name']]

library(forcats)
project_srt@meta.data[is.na(project_srt$lesions), 'lesions'] = 'Regev'
project_srt[['new_name']] = NA
project_srt@meta.data[rownames(project_srt@meta.data[project_srt$From=='CRC', ]), 'new_name'] = project_srt@meta.data[project_srt$From=='CRC', 'cluster_name']
project_srt@meta.data[rownames(project_srt@meta.data[project_srt$From=='Regev', ]), 'new_name'] = project_srt@meta.data[project_srt$From=='Regev', 'predicted.cluster_name']

project_srt=project_srt %>% subset(new_name %in% c('endo_PLVAP','iCAF_CXCL14','mCAF_RGS5'),invert=T)

data = project_srt@meta.data %>%
  filter(lesions!='Liver_meta') %>%
  group_by(orig.ident, new_name) %>%
  summarise(num=n()) 

x = reshape2::dcast(data, orig.ident~new_name)
rownames(x) = x[,1]
x = x[,-1]
x[is.na(x)] = 0
x = reshape2::melt(as.matrix(x))
colnames(x) = c('sample', 'new_name', 'cell_num')
x$Group = cluster_2_name[as.vector(x$new_name)]

x$sample_all = table(project_srt$orig.ident)[as.vector(x$sample)]

x$prop = x$cell_num / x$sample_all
x$prop[is.na(x$prop)] = 0

meta = project_srt@meta.data[, c('orig.ident', 'From')]
meta = meta[!duplicated(meta),]
rownames(meta) = meta$orig.ident

x$From = meta[as.vector(x$sample), 'From']
x = as.data.frame(x)
library(ggpubr)

FC_res =c()
unique_cell_name = unique(as.vector(x$new_name))
for(i in unique_cell_name){
  my_prop = x %>% filter(new_name==i&From=='CRC') %>% select(prop)
  he_prop = x %>% filter(new_name==i&From=='Regev') %>% select(prop)
  logFC = log2(mean(my_prop$prop) / mean(he_prop$prop))
  pvalue = wilcox.test(my_prop$prop, he_prop$prop)$p.value
  FC_res = rbind(FC_res, c(logFC, pvalue))
}
FC_res = as.data.frame(FC_res)
rownames(FC_res) = unique_cell_name
colnames(FC_res) = c('FC', 'pvalue')

drop_name = x %>% filter(From=='Regev') %>%
  group_by(new_name) %>%
  summarise(mean_prop=mean(prop)) %>% as.data.frame()

drop_name = as.vector(drop_name[drop_name$mean_prop <=0.005, 'new_name'])


FC_res[drop_name, 'pvalue'] = 1 # is.infinite(FC_res$FC)
FC_res[drop_name, 'FC'] = 1


library(ComplexHeatmap)
cell_mata = project_srt@meta.data %>%
  group_by(cluster_name) %>%
  summarise(num=n()) %>% as.data.frame()


top_anno_df = na.omit(cell_mata)
top_anno_df$Cluster='Tcell'
top_anno_df$Cluster[top_anno_df$cluster_name%in%c("Plasma_c1_IGHA2","B_c2_CXCR4")]='Bcell'
top_anno_df$Cluster[top_anno_df$cluster_name%in%c('Macro_c1_CCL18','Mono_c2_FCN1','Macro_c3_CD163L1','cDC2_c4_CLEC10A','Mono_c5_ANGPTL4',
                                                  'Cycling_Myeloid_c6','Macro_c7_CXCL10','DC_c8_LAMP3','pDC_c9_LILRA4','Macro_c10_MT1G','Macro_c11_HSPH1','Macro_c12_CD5L','cDC1_c13_CLEC9A',
                                                  'Neutrophil','Mast_cell')]='Myeloid'


rownames(top_anno_df) = top_anno_df$cluster_name
top_anno_df$FC = FC_res[rownames(top_anno_df), 'FC']
top_anno_df$pvalue = FC_res[rownames(top_anno_df), 'pvalue']

top_anno_df$pvalue2 = log(top_anno_df$pvalue )* sign(top_anno_df$FC)

top_anno_df = top_anno_df %>% arrange(Cluster, pvalue2) %>% as.data.frame()
top_anno_df$Cluster = factor(top_anno_df$Cluster, levels = c('Tcell', 'Bcell', 'Myeloid'))
top_anno_df2 = rbind(top_anno_df[top_anno_df$Cluster=='Tcell', ], top_anno_df[top_anno_df$Cluster=='Bcell', ])
top_anno_df2 = rbind(top_anno_df2, top_anno_df[top_anno_df$Cluster=='Myeloid', ])
top_anno_df = as.data.frame(top_anno_df2)

FC_res = FC_res[top_anno_df$cluster_name, ]
ht_data = t(FC_res)[1,,drop=F]
ht_data_p =  -log10(t(FC_res)[2,,drop=F]) 


col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", 'gray',  "red"))

Cluster_color=c('Tcell'='#b56cc1','Bcell'='#49a87e',
                'Myeloid'='#fec836', 'Stromal'='#ef5f5f')

p_value_col_fun <- circlize::colorRamp2(c(0, 1), c("black", "white"))
p_anno <- HeatmapAnnotation(p_value = anno_barplot(ht_data_p[1,], border = TRUE, fill = p_value_col_fun))

ht_opt(heatmap_border = TRUE)
width_in =  360/ 25.4
height_in = 100 / 25.4
dev.off()
pdf('fig 1E.pdf', width=width_in, height=height_in)

Heatmap(ht_data,
        name = "log2FC",
        cluster_rows = F, cluster_columns = F,
        col = col_fun, rect_gp = gpar(type = "none"), 
        width = ncol(ht_data)*unit(10, "mm"), # 宽度
        height = nrow(ht_data)*unit(10, "mm"), # 高度
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, 
                    height = height, 
                    gp = gpar(col = "grey", fill = NA))
          
          grid.circle(x = x, y = y, r = max(ht_data_p[i, j] * min(unit.c(width, height)), 2*min(unit.c(width, height))), 
                      gp = gpar(fill = col_fun(ht_data[i, j]), col = NA))
        },
        top_annotation = HeatmapAnnotation(df=top_anno_df[,'Cluster',drop=F],
                                           col  = list(Cluster=Cluster_color)),
        bottom_annotation = p_anno
)

dev.off()





###fig S1F####
new_level=colnames(ht_data) %>% setdiff(drop_name)

df_box=x %>% filter(!new_name%in%c(drop_name))
df_box$group[df_box$From=='CRC']='mCRC'
df_box$group[df_box$From=='Regev']='nmCRC'
df_box$group=factor(df_box$group,levels=c('nmCRC','mCRC'))

df_box$new_name=factor(df_box$new_name,levels = new_level)

df_box%>% 
  ggboxplot(x='new_name', y='prop', color='group', add='jitter',size = 0.5,shape = 16, add.params = list(size = 1),outlier.shape = NA)+
  theme(axis.text.x = element_text(angle=60, hjust=1, vjust=1, size = 10))+
  stat_compare_means(aes(group=group), label='p.signif', hide.ns = T)+ # 'p.format'
  xlab('')+ylab('Cell proportion of all immune cells')+
  scale_color_d3()


#remove outlier
library(ggplot2)
library(ggpubr)
library(dplyr)

box_stats <- df_box %>%
  group_by(new_name, group) %>%
  summarize(
    Q1 = quantile(prop, 0.25),
    Q3 = quantile(prop, 0.75)
  ) %>%
  mutate(
    lower = Q1 - 1.5 * (Q3 - Q1),
    upper = Q3 + 1.5 * (Q3 - Q1)
  )

df_no_outliers <- df_box %>%
  left_join(box_stats, by = c("new_name", "group")) %>%
  filter(prop >= lower & prop <= upper) %>%
  select(-c(Q1, Q3, lower, upper))

df_box %>%
  ggboxplot(x = 'new_name', y = 'prop', color = 'group', outlier.shape = NA, size = 0.5) +
  geom_jitter(data = df_no_outliers, 
              aes(x = new_name, y = prop, color = group, group = group), 
              position = position_jitterdodge(dodge.width = 0.8), 
              size = 1, shape = 16) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 10)) +
  stat_compare_means(aes(group = group), label = 'p.format', hide.ns = TRUE, label.y = 0.45) +
  #theme(text = element_text(size = 2)) +
  xlab('') + ylab('Cell proportion of all immune cells') +
  scale_color_d3()+ coord_cartesian(ylim = c(0, 0.5))

ggsave("fig S1F.pdf",width = 7,height = 4,units = "in")




#fig S2######
#install.packages("reticulate")
library(dplyr)
library(reticulate)
library(Seurat)
library(ggsci)
library(ggplot2)
library(SeuratDisk)

project_srt[['CellType']] = project_srt$predicted.cluster_name
project_srt@meta.data[is.na(project_srt$predicted.cluster_name), 'CellType'] = 'CRC_cells'

my_srt=project_srt %>% subset(From=='CRC')
new_srt=project_srt %>% subset(From=='Regev')
DimPlot(my_srt, reduction = "umap",group.by = 'cluster_name')+
  DimPlot(new_srt, reduction = "umap", group.by = 'predicted.cluster_name')

SaveH5Seurat(my_srt, filename = "adata_SYSU.h5seurat", overwrite=TRUE)
Convert("adata_SYSU.h5seurat", dest = "h5ad", overwrite = TRUE)

SaveH5Seurat(new_srt, filename = "adata_regev.h5seurat", overwrite=TRUE)
Convert("adata_regev.h5seurat", dest = "h5ad", overwrite = TRUE)


####CellTypist
regev_celltypist=read.csv("regev_celltypist.csv",header = T)


new_srt@meta.data$X=rownames(new_srt@meta.data)
new_srt@meta.data=merge(new_srt@meta.data,regev_celltypist,by='X',all = TRUE)
rownames(new_srt@meta.data)=new_srt@meta.data$X


new_srt$CellType=factor(new_srt$CellType,levels = c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                                    'CD8_CXCL13_Tex','CD8_GZMK_Tem','Cycling_T','TYROBP_NK&NKT','Plasma_c1_IGHA2','B_c2_CXCR4'
                                                    ,'Macro_c1_CCL18','Mono_c2_FCN1','Macro_c3_CD163L1','cDC2_c4_CLEC10A','Mono_c5_ANGPTL4',
                                                    'Cycling_Myeloid_c6','Macro_c7_CXCL10','DC_c8_LAMP3','pDC_c9_LILRA4',
                                                    'Neutrophil','Mast_cell'))

#Seurat mapping plot
DimPlot(new_srt, seed = 1, reduction = "umap", 
        group.by = "CellType",label = F, raster = F)+
  scale_color_manual(values = cell_color)+
  labs(title = 'nmCRC (Seurat)')
ggsave("nmCRC_mapping_Seurat.pdf",width = 8,height = 4,units = "in")


#mCRC mapping plot
my_srt=my_srt %>% subset(cluster_name%in%c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                           'CD8_CXCL13_Tex','CD8_GZMK_Tem','CD8_SLC4A10_MAIT','Cycling_T','TYROBP_NK&NKT','Plasma_c1_IGHA2','B_c2_CXCR4',
                                           'Macro_c1_CCL18','Mono_c2_FCN1','Macro_c3_CD163L1','cDC2_c4_CLEC10A','Mono_c5_ANGPTL4',
                                           'Cycling_Myeloid_c6','Macro_c7_CXCL10','DC_c8_LAMP3','pDC_c9_LILRA4','Macro_c10_MT1G','Macro_c11_HSPH1','Macro_c12_CD5L','cDC1_c13_CLEC9A',
                                           'Neutrophil','Mast_cell'))

my_srt$cluster_name=factor(my_srt$cluster_name,levels = c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                                          'CD8_CXCL13_Tex','CD8_GZMK_Tem','CD8_SLC4A10_MAIT','Cycling_T','TYROBP_NK&NKT','Plasma_c1_IGHA2','B_c2_CXCR4',
                                                          'Macro_c1_CCL18','Mono_c2_FCN1','Macro_c3_CD163L1','cDC2_c4_CLEC10A','Mono_c5_ANGPTL4',
                                                          'Cycling_Myeloid_c6','Macro_c7_CXCL10','DC_c8_LAMP3','pDC_c9_LILRA4','Macro_c10_MT1G','Macro_c11_HSPH1','Macro_c12_CD5L','cDC1_c13_CLEC9A',
                                                          'Neutrophil','Mast_cell'))
#Seurat mapping plot
DimPlot(my_srt, seed = 1, reduction = "umap", 
        group.by = "cluster_name",label = F, raster = F)+
  scale_color_manual(values = cell_color)+
  labs(title = 'mCRC (Seurat)')
ggsave("mCRC_mapping_Seurat.pdf",width = 8,height = 4,units = "in")




#CellTypist mapping plot
new_srt$majority_voting=factor(new_srt$majority_voting,levels = c('CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                                                  'CD8_CXCL13_Tex','CD8_GZMK_Tem','CD8_SLC4A10_MAIT','Cycling_T','TYROBP_NK&NKT','Plasma_c1_IGHA2','B_c2_CXCR4',
                                                                  'Macro_c1_CCL18','Mono_c2_FCN1','Macro_c3_CD163L1','cDC2_c4_CLEC10A','Mono_c5_ANGPTL4',
                                                                  'Cycling_Myeloid_c6','Macro_c7_CXCL10','DC_c8_LAMP3','pDC_c9_LILRA4','Macro_c10_MT1G','Macro_c11_HSPH1','cDC1_c13_CLEC9A',
                                                                  'Neutrophil','Mast_cell'))

table(new_srt$majority_voting)
unique(new_srt$majority_voting)
head(new_srt)
DimPlot(new_srt, seed = 1, reduction = "umap", 
        group.by = "majority_voting",label = F, raster = F)+
  scale_color_manual(values = cell_color)+
  labs(title = 'nmCRC (CellTypist)')
ggsave("nmCRC_mapping_CellTypist.pdf",width = 8,height = 4,units = "in")



####scANVI
regev_scANVI=read.csv("regev_scANVI.csv",header = T)


#new_srt@meta.data$X=rownames(new_srt@meta.data)
new_srt@meta.data=merge(new_srt@meta.data,regev_scANVI,by='X',all = TRUE)
rownames(new_srt@meta.data)=new_srt@meta.data$X
table(new_srt@meta.data$predictions_scanvi)



#scANVI mapping plot
new_srt$predictions_scanvi=factor(new_srt$predictions_scanvi,levels = c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                                                        'CD8_CXCL13_Tex','CD8_GZMK_Tem','CD8_SLC4A10_MAIT','Cycling_T','TYROBP_NK&NKT','Plasma_c1_IGHA2','B_c2_CXCR4',
                                                                        'Macro_c1_CCL18','Mono_c2_FCN1','Macro_c3_CD163L1','cDC2_c4_CLEC10A','Mono_c5_ANGPTL4',
                                                                        'Cycling_Myeloid_c6','Macro_c7_CXCL10','DC_c8_LAMP3','pDC_c9_LILRA4','Macro_c10_MT1G','Macro_c11_HSPH1','Macro_c12_CD5L','cDC1_c13_CLEC9A',
                                                                        'Neutrophil','Mast_cell'))

table(new_srt$predictions_scanvi)
unique(new_srt$predictions_scanvi)
head(new_srt)
DimPlot(new_srt, seed = 1, reduction = "umap", 
        group.by = "predictions_scanvi",label = F, raster = F)+
  scale_color_manual(values = cell_color)+
  labs(title = 'nmCRC (scANVI)')
ggsave("nmCRC_mapping_scANVI.pdf",width = 8,height = 4,units = "in")








head(new_srt)
head(new_srt@active.ident)
new_srt$majority_voting <- as.factor(new_srt$majority_voting)
new_srt@active.ident=new_srt$majority_voting
head(new_srt@active.ident)
Markers_CellTypist=FindAllMarkers(new_srt, group.by='majority_voting',only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Markers_CellTypist,"Markers_CellTypist.csv")


head(new_srt)
head(new_srt@active.ident)
new_srt$predictions_scanvi <- as.factor(new_srt$predictions_scanvi)
new_srt@active.ident=new_srt$predictions_scanvi
head(new_srt@active.ident)
Markers_scANVI=FindAllMarkers(new_srt, group.by='predictions_scanvi',only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Markers_scANVI,"Markers_scANVI.csv")


head(new_srt)
head(new_srt@active.ident)
new_srt$CellType <- as.factor(new_srt$CellType)
new_srt@active.ident=new_srt$CellType
head(new_srt@active.ident)
Markers_Seurat=FindAllMarkers(new_srt, group.by='CellType',only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Markers_Seurat,"Markers_Seurat.csv")


top_10_Seurat=Markers_Seurat %>% filter(pct.2<0.2) %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) 
top_10_CellTypist=Markers_CellTypist %>% filter(pct.2<0.2) %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) 
top_10_scANVI=Markers_scANVI %>% filter(pct.2<0.2) %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) 


new_srt <- NormalizeData(new_srt, normalization.method = "LogNormalize", scale.factor = 10000)
new_srt <- FindVariableFeatures(new_srt, selection.method = "vst", nfeatures = 2000)
new_srt <- ScaleData(new_srt, features = top_2_Seurat$gene)




library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

#####Seurat_marker_gene_expression######

gene_list=c("CD6","IL2RA","CCL20","IL7R","CXCL13","GZMK","STMN1","GNLY","IGHA2","MS4A1","CCL18",  
            "FCN1","CD163L1","CLEC10A","ANGPTL4","MKI67","CXCL10","LAMP3","LILRA4","G0S2","TPSB2")

avg_expression <- new_srt@assays$RNA@data[gene_list, ] %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = new_srt$CellType) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean))

avg_expression_long <- avg_expression %>%
  gather(gene, mean, -cluster)

exp_data <- new_srt@assays$RNA@data[gene_list, ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = new_srt$CellType) %>%
  gather(gene, exp, -cluster)

final_data <- left_join(avg_expression_long, exp_data, by = c("cluster", "gene"))
final_data=final_data %>% filter(!cluster%in%c('Macro_c11_HSPH1','cDC1_c13_CLEC9A'))
final_data <- final_data %>% filter(!is.na(cluster))
final_data$cluster <- droplevels(final_data$cluster)
final_data$cluster=factor(final_data$cluster,levels = c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                                        'CD8_CXCL13_Tex','CD8_GZMK_Tem','Cycling_T','TYROBP_NK&NKT','Plasma_c1_IGHA2','B_c2_CXCR4'
                                                        ,'Macro_c1_CCL18','Mono_c2_FCN1','Macro_c3_CD163L1','cDC2_c4_CLEC10A','Mono_c5_ANGPTL4',
                                                        'Cycling_Myeloid_c6','Macro_c7_CXCL10','DC_c8_LAMP3','pDC_c9_LILRA4',
                                                        'Neutrophil','Mast_cell'))
final_data$gene=factor(final_data$gene,levels = c("CD6","IL2RA","CCL20","IL7R","CXCL13","GZMK","STMN1","GNLY","IGHA2","MS4A1","CCL18",  
                                                  "FCN1","CD163L1","CLEC10A","ANGPTL4","MKI67","CXCL10","LAMP3","LILRA4","G0S2","TPSB2"))

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
  labs(title = 'marker gene expressions (Seurat)')+ theme(plot.title = element_text(hjust = 0.5))+xlab('')

ggsave("Regev_marker_gene_Seurat.pdf",width = 5,height = 5,units = "in")




#####mCRC_Seurat_marker_gene_expression######
gene_list=c("CD6","IL2RA","CCL20","IL7R","CXCL13","GZMK","STMN1","GNLY","IGHA2","MS4A1","CCL18",  
            "FCN1","CD163L1","CLEC10A","ANGPTL4","MKI67","CXCL10","LAMP3","LILRA4","G0S2","TPSB2")

avg_expression <- my_srt@assays$RNA@data[gene_list, ] %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = my_srt$cluster_name) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean))

avg_expression_long <- avg_expression %>%
  gather(gene, mean, -cluster)

exp_data <- my_srt@assays$RNA@data[gene_list, ] %>%
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = my_srt$cluster_name) %>%
  gather(gene, exp, -cluster)

final_data <- left_join(avg_expression_long, exp_data, by = c("cluster", "gene"))
unique(final_data$cluster)
final_data=final_data %>% subset(!(cluster%in%c('CD8_SLC4A10_MAIT','Macro_c10_MT1G','Macro_c11_HSPH1','Macro_c12_CD5L','cDC1_c13_CLEC9A',"endo_PLVAP","iCAF_CXCL14","mCAF_RGS5")))
final_data$cluster <- droplevels(final_data$cluster)
final_data$cluster=factor(final_data$cluster,levels = c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                                        'CD8_CXCL13_Tex','CD8_GZMK_Tem','Cycling_T','TYROBP_NK&NKT','Plasma_c1_IGHA2','B_c2_CXCR4',
                                                        'Macro_c1_CCL18','Mono_c2_FCN1','Macro_c3_CD163L1','cDC2_c4_CLEC10A','Mono_c5_ANGPTL4',
                                                        'Cycling_Myeloid_c6','Macro_c7_CXCL10','DC_c8_LAMP3','pDC_c9_LILRA4',
                                                        'Neutrophil','Mast_cell'))
final_data$gene=factor(final_data$gene,levels = c("CD6",'IL2RA','CCL20','IL7R','CXCL13','GZMK',
                                                  'STMN1','GNLY','IGHA2',"MS4A1",'CCL18',
                                                  'FCN1','CD163L1','CLEC10A','ANGPTL4','MKI67','CXCL10',
                                                  'LAMP3','LILRA4',
                                                  'G0S2','TPSB2'))

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
        axis.ticks = element_line(size=0.1))+
  labs(title = 'marker gene expressions (Seurat)')+ theme(plot.title = element_text(hjust = 0.5))+xlab('')

ggsave("mCRC_marker_gene_Seurat.pdf",width = 5.5,height = 5,units = "in")




#correlation heatmap######


gene_list=c("IL2RA","CCL20","IL7R","CXCL13","GZMK","STMN1","GNLY","IGHA2","MS4A1","CCL18",
            "FCN1","CD163L1","CLEC10A","ANGPTL4","MKI67","CXCL10","LAMP3","LILRA4","G0S2","TPSB2")

cluster_keep=intersect(unique(new_srt$majority_voting),unique(new_srt$CellType))


#CellTypist
avg_expression_CellTypist <- new_srt@assays$RNA@data[gene_list, ] %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = new_srt$majority_voting) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean)) %>% 
  filter(cluster%in%cluster_keep)
avg_expression_CellTypist=as.data.frame(avg_expression_CellTypist)
rownames(avg_expression_CellTypist)=avg_expression_CellTypist$cluster
avg_expression_CellTypist=avg_expression_CellTypist[,-1]

#Seurat_nmCRC
avg_expression_Seurat_nmCRC <- new_srt@assays$RNA@data[gene_list, ] %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = new_srt$CellType) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean))%>% 
  filter(cluster%in%cluster_keep)
avg_expression_Seurat_nmCRC=as.data.frame(avg_expression_Seurat_nmCRC)
rownames(avg_expression_Seurat_nmCRC)=avg_expression_Seurat_nmCRC$cluster
avg_expression_Seurat_nmCRC=avg_expression_Seurat_nmCRC[,-1]

#scANVI
avg_expression_scANVI <- new_srt@assays$RNA@data[gene_list, ] %>% 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  mutate(cluster = new_srt$predictions_scanvi) %>%
  group_by(cluster) %>%
  summarize(across(everything(), mean))%>% 
  filter(cluster%in%cluster_keep)
avg_expression_scANVI=as.data.frame(avg_expression_scANVI)
rownames(avg_expression_scANVI)=avg_expression_scANVI$cluster
avg_expression_scANVI=avg_expression_scANVI[,-1]


df2 <- data.matrix(avg_expression_Seurat_nmCRC)
df3 <- data.matrix(avg_expression_CellTypist)
df4 <- data.matrix(avg_expression_scANVI)


#Seurat Celltypist cor
cor_matrix <- matrix(NA, nrow(df2), nrow(df3))
rownames(cor_matrix) <- rownames(df2)
colnames(cor_matrix) <- rownames(df3)

for(i in 1:nrow(df2)){
  for(j in 1:nrow(df3)){
    cor_matrix[i,j] <- cor(df2[i,], df3[j,], method="pearson")[1]
  }
}

pheatmap(cor_matrix,
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90)
#pdf('correlation_heatmap_Seurat_CellTypist.pdf', width = 5.9, height = 5.1)




#Seurat scANVI cor
cor_matrix <- matrix(NA, nrow(df2), nrow(df4))
rownames(cor_matrix) <- rownames(df2)
colnames(cor_matrix) <- rownames(df4)

for(i in 1:nrow(df2)){
  for(j in 1:nrow(df4)){
    cor_matrix[i,j] <- cor(df2[i,], df4[j,], method="pearson")[1]
  }
}

pheatmap(cor_matrix,
         cluster_rows = F,
         cluster_cols = F,
         angle_col = 90)
#pdf('correlation_heatmap_Seurat_scANVI.pdf', width = 5.9, height = 5.1)






