library(monocle3)
library(ggplot2)
library(dplyr)
library(Seurat)
library(ggsci)

Tcell=readRDS("Tcell_in_ann_new_1227.rds") 
## Seurat obj for monocle3
CD4_T=Tcell %>% subset(cycling_sub%in%c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_Cycling_T'))
table(CD4_T@meta.data$cluster_name)

#####monocle3 pipeline
expression_matrix=GetAssayData(CD4_T,assay = 'RNA',slot = 'counts') 
cell_metadata=CD4_T@meta.data
gene_annotation=data.frame(gene_short_name=rownames(expression_matrix),row.names = rownames(expression_matrix))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)
plot_pc_variance_explained(cds)

## Step 2: Remove batch effects with cell alignment
#cds <- align_cds(cds, alignment_group = "batch")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds,reduction_method = c("UMAP"))
# plot_cells(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)
# plot_cells(cds)

p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cluster_name") + ggtitle('cds.umap') 
p1

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(Tcell, reduction = "umap")   
int.embed <- int.embed[rownames(cds.embed),]        
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="cluster_name") + ggtitle('int.umap')
p2
p = p1|p2
p


cds <- learn_graph(cds, learn_graph_control = list(rann.k=40),use_partition=F,close_loop=F)

cds=order_cells(cds,reduction_method = "UMAP")

scales::show_col(pal_igv()(9))


####fig 3D####
plot_cells(cds, reduction_method="UMAP", 
           alpha = 1,
           cell_size = 0.2,
           graph_label_size = 2,
           label_cell_groups = F,
           label_branch_points = F,
           label_roots = T,
           label_leaves = F,
           color_cells_by="cluster_name")+
  scale_color_manual(values=c('#5050FF','#CE3D32','#749B58','#802268'))+
  ggtitle('CD4+ T cell trajectory')

ggsave("Fig 3D.pdf",width = 7,height = 6,units = "in")


####fig 3E####
plot_cells(cds,alpha = 0.4,
           cell_size = 0.25,
           label_branch_points = F,
           label_roots = T,
           label_leaves = F,
           color_cells_by="pseudotime")+
  ggtitle('CD4+ T cell trajectory')

ggsave("Fig 3E.pdf",width = 6,height = 6,units = "in")




