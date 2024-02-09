library(dplyr)
library(ggplot2)
library(ggsignif)
library(ggsci)
library(Seurat)
library(pheatmap)

Total.sub=readRDS("Total_20220309.rds")
DimPlot(Total.sub, reduction = "umap",label = F, group.by = "Cluster")+scale_color_igv()

unique(Total.sub$cluster_name)
Total.sub$cluster_name[Total.sub$cluster_name=='Macro_c5_ANGPTL4']='Mono_c5_ANGPTL4' #rename
Total.sub$cluster_name=factor(Total.sub$cluster_name,levels = c('Tumor cell','CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_IL7R_Tn-like2',
                                                                'CD8_CXCL13_Tex','CD8_GZMK_Tem','CD8_SLC4A10_MAIT','Cycling_T','TYROBP_NK&NKT','Macro_c1_CCL18',
                                                                'Mono_c2_FCN1','Macro_c3_CD163L1','cDC2_c4_CLEC10A','Mono_c5_ANGPTL4','Cycling_Myeloid_c6',
                                                                'Macro_c7_CXCL10','DC_c8_LAMP3','pDC_c9_LILRA4','Macro_c10_MT1G','Macro_c11_HSPH1','Macro_c12_CD5L',
                                                                'cDC1_c13_CLEC9A','Neutrophil','Mast_cell','Plasma_c1_IGHA2','B_c2_CXCR4','mCAF_RGS5','iCAF_CXCL14',
                                                                'endo_PLVAP'))
#Fig 1B#####
DimPlot(Total.sub, reduction = "umap",label = F, group.by = "cluster_name")+scale_color_manual(values = pal_igv()(51)[c(1:25,30,33,34,36,39)])
#scales::show_col(pal_igv()(51))
ggsave("Fig 1B.pdf",width = 11,height = 7,units = "in")

head(Total.sub)

#Fig S1C####
DimPlot(Total.sub, reduction = "umap",label = F, group.by = "patient")+scale_color_igv()
ggsave("Fig S1C.pdf",width = 8,height = 7,units = "in")

#Fig S1D####
FeaturePlot(Total.sub, features=c('EPCAM','CD3E','CD79A','CD68','COL1A1'),reduction = "umap",cols = c("#DFDFDF","#FF0000"),ncol = 3) 
ggsave("Fig S1D.pdf",width = 8,height = 5,units = "in")


#Fig S1E#####
#Exclude mitochondrial and ribosomal genes based on a manually curated list
excludeGene<-read.table("excludeGene.txt",header=TRUE)
Total.sub<-Total.sub[!rownames(Total.sub) %in% excludeGene$Gene]  #29569 features across 160584 samples within 1 assay 
Idents(Total.sub) = Total.sub$cluster_name
markers=FindAllMarkers(Total.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#markers=read.csv("all_cell_markers.csv")

top_10_new<-markers %>% filter(pct.2<0.2) %>% group_by(cluster) %>% top_n(n=10, wt = avg_log2FC) ###Filter genes for analysis where pct.2 is less than 0.2 to enhance the specificity of the selected genes.
plotdf = AverageExpression(Total.sub, assays='RNA',features=unique(top_10_new$gene) , group.by = 'cluster_name', slot='data')[[1]] 

plotdf=scale(plotdf)
plotdf[plotdf>=3]=3
plotdf[plotdf < (-3)]= -3 

#breaks
bk <- c(seq(-3,-0.01,by=0.01),seq(0,3,by=0.01))
pheatmap(plotdf,cluster_rows = F,cluster_cols = F,
         show_colnames = T,angle_col = '90',fontsize_row=3,
         color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
         legend_breaks=seq(-3,3,1),
         gaps_col = c(1,10,25,27,29),
         gaps_row = c(10,75,165,183,202),
         breaks=bk)



