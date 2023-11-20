library(ggplot2)
library(viridis)
library(phateR)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggsci)

Tcell=readRDS("Tcell_in_ann_new_1227.rds")
CD4_T=Tcell %>% subset(cycling_sub%in%c('CD4_CCR7_Tn-like1','CD4_FOXP3_Treg','CD4_IL17A_Th17','CD4_Cycling_T'))

markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
top10.CD4<-markers %>% group_by(cluster) %>% filter(cluster%in%c(0,2,3,8)) %>% top_n(n=10, wt = avg_log2FC) 

#BiocManager::install("monocle")
library(monocle)
library(igraph) 

CD4_matrix<-as(as.matrix(GetAssayData(CD4_T,slot = "counts")), 'sparseMatrix')
dim(CD4_matrix)
CD4_matrix[1:6,1:7]
###fd,feature data
CD4_T_feature<-data.frame(gene_id=rownames(CD4_matrix),gene_short_name=rownames(CD4_matrix))
rownames(CD4_T_feature)<-rownames(CD4_matrix)
CD4_T_fd<-new("AnnotatedDataFrame", data = CD4_T_feature)
head(CD4_T_feature)
###pd,pheno data
sample_ann<-CD4_T@meta.data
rownames(sample_ann)<-colnames(CD4_matrix)
CD4_T_pd<-new("AnnotatedDataFrame", data =sample_ann)
head(sample_ann)
###cds
CD4_T.cds<-newCellDataSet(CD4_matrix,phenoData =CD4_T_pd,featureData =CD4_T_fd,lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())
head(pData(CD4_T.cds))
head(fData(CD4_T.cds))


CD4_T.cds <- estimateSizeFactors(CD4_T.cds)
CD4_T.cds <- estimateDispersions(CD4_T.cds)
CD4_T.cds<- detectGenes(CD4_T.cds, min_expr = 0.1)
CD4_T.cds<- setOrderingFilter(CD4_T.cds,ordering_genes = top10.CD4$gene)

#DDRtree
CD4_T.cds <- reduceDimension(
  CD4_T.cds,
  max_components = 2, 
  method = 'DDRTree')
CD4_T.cds <- orderCells(CD4_T.cds)

head(pData(CD4_T.cds))
CD4_T.cds <- orderCells(CD4_T.cds, root_state = 4) #reset

plot_cell_trajectory(CD4_T.cds,cell_size = 1)+scale_color_igv()


allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB",
            "#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
            "#32CD32","#F0E68C","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887","#FFFFE0")

plot_cell_trajectory(CD4_T.cds,cell_size = 1,color_by = "location")+scale_color_manual(values = allcolour)
plot_cell_trajectory(CD4_T.cds,cell_size = 1,color_by = "lesions")+scale_color_manual(values = allcolour)

plot_cell_trajectory(CD4_T.cds[CD4_T.cds$cluster_name=='CD4_IL17A_Th17'],cell_size = 1,color_by = "cluster_name")+scale_color_manual(values = allcolour)

plot_cell_trajectory(CD4_T.cds,cell_size = 1,color_by = "cluster_name")+scale_color_manual(values = c("#DC143C00","#0000FF","#20B2AA00","#FFA50000"))
plot_cell_trajectory(CD4_T.cds,cell_size = 1,color_by = "cluster_name")+scale_color_manual(values = c("#DC143C00","#0000FF00","#20B2AA","#FFA50000"))

plot_genes_in_pseudotime(CD4_T.cds)

plot_genes_branched_pseudotime(CD4_T.cds[c('IL17A','FOXP3','MKI67'),], 
                               branch_point = 1,
                               color_by = 'cluster_name',
                               ncol = 1)+scale_color_manual(values = allcolour)



####Fig 3J####
BEAM_res<-BEAM(CD4_T.cds, branch_point = 1,cores = 1)

t=plot_genes_branched_heatmap(CD4_T.cds[row.names(subset(BEAM_res, qval < 1e-4)),],
                              branch_point = 1,   
                              num_clusters = 6,   #same as tumor cells, 6 clusters
                              cores = 4,
                              use_gene_short_name = T,
                              show_rownames = F,
                              return_heatmap=T)



plot_genes_branched_heatmap(CD4_T.cds[row.names(subset(BEAM_res, qval < 1e-4)),],
                            branch_point = 1,   
                            num_clusters = 6,   
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = F)




table(t$annotation_row$Cluster)
rownames(t$annotation_row$Cluster)
branch_gene=t$annotation_row
branch_gene_name_Left=rownames(branch_gene %>% filter(Cluster=='4')) #left branch
branch_gene_name_Right=rownames(branch_gene %>% filter(Cluster=='3')) #right branch
branch_gene_name_Mid=rownames(branch_gene %>% filter(Cluster=='1')) #middle

plot_genes_branched_heatmap(CD4_T.cds[c(branch_gene_name_Left,branch_gene_name_Right),],
                            branch_point = 1,   
                            num_clusters = 2,   
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = F)


pdf("Fig 3J.pdf",width =4,height = 8)  
plot_genes_branched_heatmap(CD4_T.cds[c(branch_gene_name_Left,branch_gene_name_Mid,branch_gene_name_Right),],
                            branch_point = 1, 
                            num_clusters = 3,   
                            cores = 4,
                            use_gene_short_name = T,
                            show_rownames = F)
dev.off()

library(clusterProfiler)
library(msigdf)
H=msigdf.human %>% filter(category_code=='h') %>% select(geneset,symbol) %>% as.data.frame() #hallmark


branch_L=enricher(branch_gene_name_Left, TERM2GENE = H)
branch_R=enricher(branch_gene_name_Right, TERM2GENE = H)
branch_M=enricher(branch_gene_name_Mid, TERM2GENE = H)


head(branch)
dotplot(branch_L, showCategory=10, title='Th17->Cycling T') 
dotplot(branch_R, showCategory=10, title='Th17->Treg') 
dotplot(branch_M, showCategory=10, title='Th17') 


write.csv(branch_L@result,"enriched_genes_upper.csv")
write.csv(branch_R@result,"enriched_genes_mid.csv")
write.csv(branch_M@result,"enriched_genes_lower.csv")


