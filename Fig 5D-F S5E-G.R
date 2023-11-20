library(monocle)
library(igraph)
library(Seurat)
library(dplyr)
library(ggsci)
library(cowplot)
library(tidyverse)
library(philentropy)
library(msigdf)
library(clusterProfiler)
library(enrichplot)
library(pheatmap)



###Fig 5D and S5E####
EMT<-read.delim('EMT.txt')
EMT<-EMT[-1,]

Tumor = readRDS("tumors_filter_by_cnv_man.rds")

IG_name<-rownames(Tumor)[grep("^IG[HLK]",rownames(Tumor))]
HSP_name<-rownames(Tumor)[grep("^HSP",rownames(Tumor))]
keep_name<-setdiff(rownames(Tumor), c(IG_name,"JCHAIN")) %>% setdiff(HSP_name)
Tumor<-subset(Tumor, features=keep_name)

excludeGene<-read.table("excludeGene.txt",header=TRUE)
Tumor<-Tumor[!rownames(Tumor) %in% excludeGene$Gene]


set.seed(1)

EMT.list<-list(EMT)
Tumor<-AddModuleScore(object = Tumor,features = EMT.list, name='EMT')
Tumor$EMT_score<-Tumor$EMT1

pseudotime<-function (cds_subset, cluster_rows = TRUE, hclust_method = "ward.D2", 
                      num_clusters = 6, hmcols = NULL, add_annotation_row = NULL, 
                      add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE, 
                      norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3, 
                      trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE, 
                      cores = 1) 
{
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                         max(pData(cds_subset)$Pseudotime), length.out = 100))
  m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula, 
                       relative_expr = T, new_data = newdata)
  m = m[!apply(m, 1, sum) == 0, ]
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) == 
      FALSE) {
    m = vstExprs(cds_subset, expr_matrix = m)
  }
  else if (norm_method == "log") {
    m = log10(m + pseudocount)
  }
  m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m
  
  return(heatmap_matrix)
}



####Pt02######
Pt02<-Tumor %>% subset(patient%in%'Pt02')
Pt02_matrix<-as(as.matrix(GetAssayData(Pt02,slot = "counts")), 'sparseMatrix')
Pt02_feature<-data.frame(gene_id=rownames(Pt02_matrix),gene_short_name=rownames(Pt02_matrix))
rownames(Pt02_feature)<-rownames(Pt02_matrix)
Pt02_fd<-new("AnnotatedDataFrame", data = Pt02_feature)
sample_ann<-Pt02@meta.data
rownames(sample_ann)<-colnames(Pt02_matrix)
Pt02_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt02.cds<-newCellDataSet(Pt02_matrix,phenoData =Pt02_pd,featureData =Pt02_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt02.cds <- estimateSizeFactors(Pt02.cds)
Pt02.cds <- estimateDispersions(Pt02.cds)
Pt02.cds<- detectGenes(Pt02.cds, min_expr = 0.1)
Pt02.cds<- setOrderingFilter(Pt02.cds,ordering_genes = EMT)
Pt02.cds <- reduceDimension(
  Pt02.cds,
  max_components = 2,
  method = 'DDRTree')
Pt02.cds <- orderCells(Pt02.cds)
plot_cell_trajectory(Pt02.cds,cell_size = 0.5,color_by = "Pseudotime")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(3,'mm'))+labs(title = 'EMT')+coord_fixed(ratio = 2)+NoLegend() #+ylim(c(-1,2))

dim(Pt02.cds@reducedDimS)
plot(Pt02.cds@reducedDimS[1,],Pt02.cds@reducedDimS[2,])
plot(Pt02.cds@reducedDimW[1,],Pt02.cds@reducedDimW[2,])
sort(Pt02.cds@reducedDimS[2,]) %>% tail(10)
names(sort(Pt02.cds@reducedDimS[2,]) %>% tail(10))

Pt02.filter<-subset(Pt02,cells=names(sort(Pt02.cds@reducedDimS[2,]) %>% tail(10)),invert=T)

Pt02_matrix<-as(as.matrix(GetAssayData(Pt02.filter,slot = "counts")), 'sparseMatrix')

Pt02_feature<-data.frame(gene_id=rownames(Pt02_matrix),gene_short_name=rownames(Pt02_matrix))
rownames(Pt02_feature)<-rownames(Pt02_matrix)
Pt02_fd<-new("AnnotatedDataFrame", data = Pt02_feature)
sample_ann<-Pt02.filter@meta.data
rownames(sample_ann)<-colnames(Pt02_matrix)
Pt02_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt02.cds<-newCellDataSet(Pt02_matrix,phenoData =Pt02_pd,featureData =Pt02_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt02.cds <- estimateSizeFactors(Pt02.cds)
Pt02.cds <- estimateDispersions(Pt02.cds)
Pt02.cds<- detectGenes(Pt02.cds, min_expr = 0.1)
Pt02.cds<- setOrderingFilter(Pt02.cds,ordering_genes = EMT)
Pt02.cds <- reduceDimension(
  Pt02.cds,
  max_components = 2, 
  method = 'DDRTree')
Pt02.cds <- orderCells(Pt02.cds)

pData(Pt02.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt02.cds <- orderCells(Pt02.cds, root_state = 5)
pData(Pt02.cds)$PseudoEMT<-pData(Pt02.cds)$Pseudotime

plot_cell_trajectory(Pt02.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt02_EMTscore.pdf",width = 8,height = 7,units = "in")


plot_cell_trajectory(Pt02.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt02_PseudoEMT.pdf",width = 10,height = 8,units = "in")

diff_test_res.Pt02.total <- differentialGeneTest(Pt02.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt02.filter<-subset(diff_test_res.Pt02.total, qval < 0.05)
sig_gene_names.Pt02 <- row.names(subset(diff_test_res.Pt02.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt02.cds[sig_gene_names.Pt02,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt02.cds[sig_gene_names.Pt02,],  
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt02.filter=left_join(tep, diff_test_res.Pt02.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt02.filter$Cluster<-factor(diff_test_res.Pt02.filter$Cluster,levels = c(5,3,6,4,2,1), ordered=T) 
b<-diff_test_res.Pt02.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt02.filter$Cluster))), 
         main='Pt02 PseudoEMT DEGs',
         filename="Pt02 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()



msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(5,3,6,4,2,1)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt02=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt02)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt02=diff_test_res.Pt02.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt02$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt02'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result) 
  rownames(res@result)=NULL
  re.Pt02=rbind(re.Pt02, res@result)
}


Pt02.order=data.frame(Cluster=c(5,3,6,4,2,1), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt02.filter$Cluster=as.numeric(diff_test_res.Pt02.filter$Cluster)
diff_test_res.Pt02.filter=diff_test_res.Pt02.filter %>% left_join(Pt02.order,by='Cluster')
diff_test_res.Pt02.filter$Pt='Pt02'


####Pt01######
Pt01<-Tumor %>% subset(patient%in%'Pt01')
Pt01_matrix<-as(as.matrix(GetAssayData(Pt01,slot = "counts")), 'sparseMatrix')
Pt01_feature<-data.frame(gene_id=rownames(Pt01_matrix),gene_short_name=rownames(Pt01_matrix))
rownames(Pt01_feature)<-rownames(Pt01_matrix)
Pt01_fd<-new("AnnotatedDataFrame", data = Pt01_feature)
sample_ann<-Pt01@meta.data
rownames(sample_ann)<-colnames(Pt01_matrix)
Pt01_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt01.cds<-newCellDataSet(Pt01_matrix,phenoData =Pt01_pd,featureData =Pt01_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt01.cds <- estimateSizeFactors(Pt01.cds)
Pt01.cds <- estimateDispersions(Pt01.cds)
Pt01.cds<- detectGenes(Pt01.cds, min_expr = 0.1)
Pt01.cds<- setOrderingFilter(Pt01.cds,ordering_genes = EMT)
Pt01.cds <- reduceDimension(
  Pt01.cds,
  max_components = 2,
  method = 'DDRTree')
Pt01.cds <- orderCells(Pt01.cds)

pData(Pt01.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt01.cds <- orderCells(Pt01.cds, root_state = 2)
pData(Pt01.cds)$PseudoEMT<-pData(Pt01.cds)$Pseudotime


plot_cell_trajectory(Pt01.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt01_EMTscore.pdf",width = 8,height = 7,units = "in")


plot_cell_trajectory(Pt01.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt01_PseudoEMT.pdf",width = 10,height = 8,units = "in")


diff_test_res.Pt01.total <- differentialGeneTest(Pt01.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt01.filter<-subset(diff_test_res.Pt01.total, qval < 0.05)
sig_gene_names.Pt01 <- row.names(subset(diff_test_res.Pt01.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt01.cds[sig_gene_names.Pt01,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt01.cds[sig_gene_names.Pt01,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt01.filter=left_join(tep, diff_test_res.Pt01.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt01.filter$Cluster<-factor(diff_test_res.Pt01.filter$Cluster,levels = c(1,4,3,6,5,2), ordered=T) 
b<-diff_test_res.Pt01.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt01.filter$Cluster))), 
         main='Pt01 PseudoEMT DEGs',
         filename="Pt01 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()





msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(1,4,3,6,5,2)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt01=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt01)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt01=diff_test_res.Pt01.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt01$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt01'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt01=rbind(re.Pt01, res@result)
}

Pt01.order=data.frame(Cluster=c(1,4,3,6,5,2), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt01.filter$Cluster=as.numeric(diff_test_res.Pt01.filter$Cluster)
diff_test_res.Pt01.filter=diff_test_res.Pt01.filter %>% left_join(Pt01.order,by='Cluster')
diff_test_res.Pt01.filter$Pt='Pt01'


####Pt03######
Pt03<-Tumor %>% subset(patient%in%'Pt03')
Pt03_matrix<-as(as.matrix(GetAssayData(Pt03,slot = "counts")), 'sparseMatrix')
Pt03_feature<-data.frame(gene_id=rownames(Pt03_matrix),gene_short_name=rownames(Pt03_matrix))
rownames(Pt03_feature)<-rownames(Pt03_matrix)
Pt03_fd<-new("AnnotatedDataFrame", data = Pt03_feature)
sample_ann<-Pt03@meta.data
rownames(sample_ann)<-colnames(Pt03_matrix)
Pt03_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt03.cds<-newCellDataSet(Pt03_matrix,phenoData =Pt03_pd,featureData =Pt03_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt03.cds <- estimateSizeFactors(Pt03.cds)
Pt03.cds <- estimateDispersions(Pt03.cds)
Pt03.cds<- detectGenes(Pt03.cds, min_expr = 0.1)
Pt03.cds<- setOrderingFilter(Pt03.cds,ordering_genes = EMT)
Pt03.cds <- reduceDimension(
  Pt03.cds,
  max_components = 2,
  method = 'DDRTree')
Pt03.cds <- orderCells(Pt03.cds)
plot_cell_trajectory(Pt03.cds,cell_size = 0.5,color_by = "Pseudotime")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(3,'mm'))+labs(title = 'EMT')+coord_fixed(ratio = 2)+NoLegend() #+ylim(c(-1,2))

dim(Pt03.cds@reducedDimS)
plot(Pt03.cds@reducedDimS[1,],Pt03.cds@reducedDimS[2,])
plot(Pt03.cds@reducedDimW[1,],Pt03.cds@reducedDimW[2,])
sort(Pt03.cds@reducedDimS[2,]) %>% tail(10)
names(sort(Pt03.cds@reducedDimS[2,]) %>% tail(10))
names(sort(Pt03.cds@reducedDimS[2,]) %>% tail(3))

Pt03.filter<-subset(Pt03,cells=names(sort(Pt03.cds@reducedDimS[2,]) %>% tail(3)),invert=T)

Pt03_matrix<-as(as.matrix(GetAssayData(Pt03.filter,slot = "counts")), 'sparseMatrix')
Pt03_feature<-data.frame(gene_id=rownames(Pt03_matrix),gene_short_name=rownames(Pt03_matrix))
rownames(Pt03_feature)<-rownames(Pt03_matrix)
Pt03_fd<-new("AnnotatedDataFrame", data = Pt03_feature)
sample_ann<-Pt03.filter@meta.data
rownames(sample_ann)<-colnames(Pt03_matrix)
Pt03_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt03.cds<-newCellDataSet(Pt03_matrix,phenoData =Pt03_pd,featureData =Pt03_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt03.cds <- estimateSizeFactors(Pt03.cds)
Pt03.cds <- estimateDispersions(Pt03.cds)
Pt03.cds<- detectGenes(Pt03.cds, min_expr = 0.1)
Pt03.cds<- setOrderingFilter(Pt03.cds,ordering_genes = EMT)
Pt03.cds <- reduceDimension(
  Pt03.cds,
  max_components = 2, 
  method = 'DDRTree')

Pt03.cds <- orderCells(Pt03.cds)

head(pData(Pt03.cds))
pData(Pt03.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt03.cds <- orderCells(Pt03.cds, root_state = 7)
pData(Pt03.cds)$PseudoEMT<-pData(Pt03.cds)$Pseudotime

plot_cell_trajectory(Pt03.cds,cell_size = 0.5)+scale_color_igv()+coord_fixed(ratio = 1)
plot_cell_trajectory(Pt03.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt03_EMTscore.pdf",width = 8,height = 7,units = "in")


plot_cell_trajectory(Pt03.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt03_PseudoEMT.pdf",width = 10,height = 8,units = "in")


diff_test_res.Pt03.total <- differentialGeneTest(Pt03.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt03.filter<-subset(diff_test_res.Pt03.total, qval < 0.05)
sig_gene_names.Pt03 <- row.names(subset(diff_test_res.Pt03.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt03.cds[sig_gene_names.Pt03,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt03.cds[sig_gene_names.Pt03,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt03.filter=left_join(tep, diff_test_res.Pt03.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt03.filter$Cluster<-factor(diff_test_res.Pt03.filter$Cluster,levels = c(2,1,5,4,6,3), ordered=T) 
b<-diff_test_res.Pt03.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt03.filter$Cluster))), 
         main='Pt03 PseudoEMT DEGs',
         filename="Pt03 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()



msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(2,1,5,4,6,3)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt03=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt03)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt03=diff_test_res.Pt03.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt03$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt03'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt03=rbind(re.Pt03, res@result)
}

Pt03.order=data.frame(Cluster=c(2,1,5,4,6,3), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt03.filter$Cluster=as.numeric(diff_test_res.Pt03.filter$Cluster)
diff_test_res.Pt03.filter=diff_test_res.Pt03.filter %>% left_join(Pt03.order,by='Cluster')
diff_test_res.Pt03.filter$Pt='Pt03'


####Pt04######
Pt04<-Tumor %>% subset(patient%in%'Pt04')
Pt04_matrix<-as(as.matrix(GetAssayData(Pt04,slot = "counts")), 'sparseMatrix')
Pt04_feature<-data.frame(gene_id=rownames(Pt04_matrix),gene_short_name=rownames(Pt04_matrix))
rownames(Pt04_feature)<-rownames(Pt04_matrix)
Pt04_fd<-new("AnnotatedDataFrame", data = Pt04_feature)
sample_ann<-Pt04@meta.data
rownames(sample_ann)<-colnames(Pt04_matrix)
Pt04_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt04.cds<-newCellDataSet(Pt04_matrix,phenoData =Pt04_pd,featureData =Pt04_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt04.cds <- estimateSizeFactors(Pt04.cds)
Pt04.cds <- estimateDispersions(Pt04.cds)
Pt04.cds<- detectGenes(Pt04.cds, min_expr = 0.1)
Pt04.cds<- setOrderingFilter(Pt04.cds,ordering_genes = EMT)
Pt04.cds <- reduceDimension(
  Pt04.cds,
  max_components = 2,
  method = 'DDRTree')
Pt04.cds <- orderCells(Pt04.cds)
plot_cell_trajectory(Pt04.cds,cell_size = 0.5,color_by = "Pseudotime")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(3,'mm'))+labs(title = 'EMT')+coord_fixed(ratio = 2)+NoLegend() #+ylim(c(-1,2))

dim(Pt04.cds@reducedDimS)
plot(Pt04.cds@reducedDimS[1,],Pt04.cds@reducedDimS[2,])
plot(Pt04.cds@reducedDimW[1,],Pt04.cds@reducedDimW[2,])
sort(Pt04.cds@reducedDimS[2,]) %>% head(10)
names(sort(Pt04.cds@reducedDimS[2,]) %>% head(10))
names(sort(Pt04.cds@reducedDimS[2,]) %>% head(3))

Pt04.filter<-subset(Pt04,cells=names(sort(Pt04.cds@reducedDimS[2,]) %>% head(3)),invert=T)

Pt04_matrix<-as(as.matrix(GetAssayData(Pt04.filter,slot = "counts")), 'sparseMatrix')
Pt04_feature<-data.frame(gene_id=rownames(Pt04_matrix),gene_short_name=rownames(Pt04_matrix))
rownames(Pt04_feature)<-rownames(Pt04_matrix)
Pt04_fd<-new("AnnotatedDataFrame", data = Pt04_feature)
sample_ann<-Pt04.filter@meta.data
rownames(sample_ann)<-colnames(Pt04_matrix)
Pt04_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt04.cds<-newCellDataSet(Pt04_matrix,phenoData =Pt04_pd,featureData =Pt04_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt04.cds <- estimateSizeFactors(Pt04.cds)
Pt04.cds <- estimateDispersions(Pt04.cds)
Pt04.cds<- detectGenes(Pt04.cds, min_expr = 0.1)
Pt04.cds<- setOrderingFilter(Pt04.cds,ordering_genes = EMT)
Pt04.cds <- reduceDimension(
  Pt04.cds,
  max_components = 2,
  method = 'DDRTree')

Pt04.cds <- orderCells(Pt04.cds)

pData(Pt04.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt04.cds <- orderCells(Pt04.cds, root_state = 6)
pData(Pt04.cds)$PseudoEMT<-pData(Pt04.cds)$Pseudotime

plot_cell_trajectory(Pt04.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt04_EMTscore.pdf",width = 8,height = 7,units = "in")

plot_cell_trajectory(Pt04.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt04_PseudoEMT.pdf",width = 10,height = 8,units = "in")

diff_test_res.Pt04.total <- differentialGeneTest(Pt04.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt04.filter<-subset(diff_test_res.Pt04.total, qval < 0.05)
sig_gene_names.Pt04 <- row.names(subset(diff_test_res.Pt04.total, qval < 0.05)) 

t <- plot_pseudotime_heatmap(Pt04.cds[sig_gene_names.Pt04,],
                             num_clusters = 6, 
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt04.cds[sig_gene_names.Pt04,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt04.filter=left_join(tep, diff_test_res.Pt04.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt04.filter$Cluster<-factor(diff_test_res.Pt04.filter$Cluster,levels = c(4,6,5,3,2,1), ordered=T) 
b<-diff_test_res.Pt04.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt04.filter$Cluster))), 
         main='Pt04 PseudoEMT DEGs',
         filename="Pt04 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()



msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(4,6,5,3,2,1)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt04=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt04)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt04=diff_test_res.Pt04.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt04$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt04'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt04=rbind(re.Pt04, res@result)
}

Pt04.order=data.frame(Cluster=c(4,6,5,3,2,1), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt04.filter$Cluster=as.numeric(diff_test_res.Pt04.filter$Cluster)
diff_test_res.Pt04.filter=diff_test_res.Pt04.filter %>% left_join(Pt04.order,by='Cluster')
diff_test_res.Pt04.filter$Pt='Pt04'


####Pt18######
Pt18<-Tumor %>% subset(patient%in%'Pt18')
Pt18_matrix<-as(as.matrix(GetAssayData(Pt18,slot = "counts")), 'sparseMatrix')
Pt18_feature<-data.frame(gene_id=rownames(Pt18_matrix),gene_short_name=rownames(Pt18_matrix))
rownames(Pt18_feature)<-rownames(Pt18_matrix)
Pt18_fd<-new("AnnotatedDataFrame", data = Pt18_feature)
sample_ann<-Pt18@meta.data
rownames(sample_ann)<-colnames(Pt18_matrix)
Pt18_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt18.cds<-newCellDataSet(Pt18_matrix,phenoData =Pt18_pd,featureData =Pt18_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt18.cds <- estimateSizeFactors(Pt18.cds)
Pt18.cds <- estimateDispersions(Pt18.cds)
Pt18.cds<- detectGenes(Pt18.cds, min_expr = 0.1)
Pt18.cds<- setOrderingFilter(Pt18.cds,ordering_genes = EMT)
Pt18.cds <- reduceDimension(
  Pt18.cds,
  max_components = 2,
  method = 'DDRTree')
Pt18.cds <- orderCells(Pt18.cds)

head(pData(Pt18.cds))
pData(Pt18.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt18.cds <- orderCells(Pt18.cds, root_state = 1)
pData(Pt18.cds)$PseudoEMT<-pData(Pt18.cds)$Pseudotime

plot_cell_trajectory(Pt18.cds,cell_size = 0.5)+scale_color_igv()+coord_fixed(ratio = 1)
plot_cell_trajectory(Pt18.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt18_EMTscore.pdf",width = 8,height = 7,units = "in")


plot_cell_trajectory(Pt18.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt18_PseudoEMT.pdf",width = 10,height = 8,units = "in")


diff_test_res.Pt18.total <- differentialGeneTest(Pt18.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt18.filter<-subset(diff_test_res.Pt18.total, qval < 0.05)
sig_gene_names.Pt18 <- row.names(subset(diff_test_res.Pt18.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt18.cds[sig_gene_names.Pt18,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt18.cds[sig_gene_names.Pt18,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt18.filter=left_join(tep, diff_test_res.Pt18.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt18.filter$Cluster<-factor(diff_test_res.Pt18.filter$Cluster,levels = c(3,1,4,2,5,6), ordered=T) 
b<-diff_test_res.Pt18.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt18.filter$Cluster))), 
         main='Pt18 PseudoEMT DEGs',
         filename="Pt18 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()



msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(3,1,4,2,5,6)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt18=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt18)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt18=diff_test_res.Pt18.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt18$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt18'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt18=rbind(re.Pt18, res@result)
}

Pt18.order=data.frame(Cluster=c(3,1,4,2,5,6), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt18.filter$Cluster=as.numeric(diff_test_res.Pt18.filter$Cluster)
diff_test_res.Pt18.filter=diff_test_res.Pt18.filter %>% left_join(Pt18.order,by='Cluster')
diff_test_res.Pt18.filter$Pt='Pt18'



#merge Pt01-04, Pt18#####
Hallmark.part1=re.Pt01 %>% rbind(re.Pt02) %>% rbind(re.Pt03) %>% rbind(re.Pt04) %>% rbind(re.Pt18)
#saveRDS(Hallmark.part1,"Hallmark_part1_1127.rds")

diff_test_res.part1=diff_test_res.Pt01.filter %>% 
  rbind(diff_test_res.Pt02.filter) %>% 
  rbind(diff_test_res.Pt03.filter) %>% 
  rbind(diff_test_res.Pt04.filter) %>% 
  rbind(diff_test_res.Pt18.filter)
#saveRDS(diff_test_res.part1,"diff_test_res_part1_1127.rds")


pData.part1=pData(Pt01.cds) %>% rbind(pData(Pt02.cds)) %>% rbind(pData(Pt03.cds)) %>% rbind(pData(Pt04.cds)) %>% rbind(pData(Pt18.cds))
#saveRDS(pData.part1, "pData_part1_1129.rds")





####Pt05######
Pt05<-Tumor %>% subset(patient%in%'Pt05')
Pt05_matrix<-as(as.matrix(GetAssayData(Pt05,slot = "counts")), 'sparseMatrix')
Pt05_feature<-data.frame(gene_id=rownames(Pt05_matrix),gene_short_name=rownames(Pt05_matrix))
rownames(Pt05_feature)<-rownames(Pt05_matrix)
Pt05_fd<-new("AnnotatedDataFrame", data = Pt05_feature)
sample_ann<-Pt05@meta.data
rownames(sample_ann)<-colnames(Pt05_matrix)
Pt05_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt05.cds<-newCellDataSet(Pt05_matrix,phenoData =Pt05_pd,featureData =Pt05_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt05.cds <- estimateSizeFactors(Pt05.cds)
Pt05.cds <- estimateDispersions(Pt05.cds)
Pt05.cds<- detectGenes(Pt05.cds, min_expr = 0.1)
Pt05.cds<- setOrderingFilter(Pt05.cds,ordering_genes = EMT)
Pt05.cds <- reduceDimension(
  Pt05.cds,
  max_components = 2, 
  method = 'DDRTree')
Pt05.cds <- orderCells(Pt05.cds)

pData(Pt05.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt05.cds <- orderCells(Pt05.cds, root_state = 2)
pData(Pt05.cds)$PseudoEMT<-pData(Pt05.cds)$Pseudotime

plot_cell_trajectory(Pt05.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt05_EMTscore.pdf",width = 8,height = 7,units = "in")

plot_cell_trajectory(Pt05.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt05_PseudoEMT.pdf",width = 10,height = 8,units = "in")


diff_test_res.Pt05.total <- differentialGeneTest(Pt05.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt05.filter<-subset(diff_test_res.Pt05.total, qval < 0.05)
sig_gene_names.Pt05 <- row.names(subset(diff_test_res.Pt05.total, qval < 0.05)) 

t <- plot_pseudotime_heatmap(Pt05.cds[sig_gene_names.Pt05,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt05.cds[sig_gene_names.Pt05,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt05.filter=left_join(tep, diff_test_res.Pt05.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt05.filter$Cluster<-factor(diff_test_res.Pt05.filter$Cluster,levels = c(4,1,2,5,6,3), ordered=T) 
b<-diff_test_res.Pt05.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt05.filter$Cluster))), 
         main='Pt05 PseudoEMT DEGs',
         filename="Pt05 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()



msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(4,1,2,5,6,3)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt05=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt05)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt05=diff_test_res.Pt05.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt05$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt05'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt05=rbind(re.Pt05, res@result)
}

Pt05.order=data.frame(Cluster=c(4,1,2,5,6,3), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt05.filter$Cluster=as.numeric(diff_test_res.Pt05.filter$Cluster)
diff_test_res.Pt05.filter=diff_test_res.Pt05.filter %>% left_join(Pt05.order,by='Cluster')
diff_test_res.Pt05.filter$Pt='Pt05'


####Pt06######
Pt06<-Tumor %>% subset(patient%in%'Pt06')
Pt06_matrix<-as(as.matrix(GetAssayData(Pt06,slot = "counts")), 'sparseMatrix')
Pt06_feature<-data.frame(gene_id=rownames(Pt06_matrix),gene_short_name=rownames(Pt06_matrix))
rownames(Pt06_feature)<-rownames(Pt06_matrix)
Pt06_fd<-new("AnnotatedDataFrame", data = Pt06_feature)
sample_ann<-Pt06@meta.data
rownames(sample_ann)<-colnames(Pt06_matrix)
Pt06_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt06.cds<-newCellDataSet(Pt06_matrix,phenoData =Pt06_pd,featureData =Pt06_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt06.cds <- estimateSizeFactors(Pt06.cds)
Pt06.cds <- estimateDispersions(Pt06.cds)
Pt06.cds<- detectGenes(Pt06.cds, min_expr = 0.1)
Pt06.cds<- setOrderingFilter(Pt06.cds,ordering_genes = EMT)
Pt06.cds <- reduceDimension(
  Pt06.cds,
  max_components = 2, 
  method = 'DDRTree')
Pt06.cds <- orderCells(Pt06.cds)
pData(Pt06.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt06.cds <- orderCells(Pt06.cds, root_state = 2)
pData(Pt06.cds)$PseudoEMT<-pData(Pt06.cds)$Pseudotime

plot_cell_trajectory(Pt06.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt06_EMTscore.pdf",width = 8,height = 7,units = "in")

plot_cell_trajectory(Pt06.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt06_PseudoEMT.pdf",width = 10,height = 8,units = "in")


diff_test_res.Pt06.total <- differentialGeneTest(Pt06.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt06.filter<-subset(diff_test_res.Pt06.total, qval < 0.05)
sig_gene_names.Pt06 <- row.names(subset(diff_test_res.Pt06.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt06.cds[sig_gene_names.Pt06,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt06.cds[sig_gene_names.Pt06,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt06.filter=left_join(tep, diff_test_res.Pt06.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt06.filter$Cluster<-factor(diff_test_res.Pt06.filter$Cluster,levels = c(5,3,4,6,1,2), ordered=T) 
b<-diff_test_res.Pt06.filter %>% arrange(Cluster)

pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt06.filter$Cluster))), 
         main='Pt06 PseudoEMT DEGs',
         filename="Pt06 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()




msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(5,3,4,6,1,2)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt06=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt06)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt06=diff_test_res.Pt06.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt06$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_')
  res@result$Pt='Pt06'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt06=rbind(re.Pt06, res@result)
}


Pt06.order=data.frame(Cluster=c(5,3,4,6,1,2), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt06.filter$Cluster=as.numeric(diff_test_res.Pt06.filter$Cluster)
diff_test_res.Pt06.filter=diff_test_res.Pt06.filter %>% left_join(Pt06.order,by='Cluster')
diff_test_res.Pt06.filter$Pt='Pt06'



####Pt07######
Pt07<-Tumor %>% subset(patient%in%'Pt07')
Pt07_matrix<-as(as.matrix(GetAssayData(Pt07,slot = "counts")), 'sparseMatrix')
Pt07_feature<-data.frame(gene_id=rownames(Pt07_matrix),gene_short_name=rownames(Pt07_matrix))
rownames(Pt07_feature)<-rownames(Pt07_matrix)
Pt07_fd<-new("AnnotatedDataFrame", data = Pt07_feature)
sample_ann<-Pt07@meta.data
rownames(sample_ann)<-colnames(Pt07_matrix)
Pt07_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt07.cds<-newCellDataSet(Pt07_matrix,phenoData =Pt07_pd,featureData =Pt07_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt07.cds <- estimateSizeFactors(Pt07.cds)
Pt07.cds <- estimateDispersions(Pt07.cds)
Pt07.cds<- detectGenes(Pt07.cds, min_expr = 0.1)
Pt07.cds<- setOrderingFilter(Pt07.cds,ordering_genes = EMT)
Pt07.cds <- reduceDimension(
  Pt07.cds,
  max_components = 2,
  method = 'DDRTree')
Pt07.cds <- orderCells(Pt07.cds)
plot_cell_trajectory(Pt07.cds,cell_size = 0.5,color_by = "Pseudotime")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(3,'mm'))+labs(title = 'EMT')+coord_fixed(ratio = 2)+NoLegend() #+ylim(c(-1,2))

dim(Pt07.cds@reducedDimS)
plot(Pt07.cds@reducedDimS[1,],Pt07.cds@reducedDimS[2,])
plot(Pt07.cds@reducedDimW[1,],Pt07.cds@reducedDimW[2,])
sort(Pt07.cds@reducedDimS[2,]) %>% tail(10)
names(sort(Pt07.cds@reducedDimS[2,]) %>% tail(2))

Pt07.filter<-subset(Pt07,cells=names(sort(Pt07.cds@reducedDimS[2,]) %>% tail(2)),invert=T)

Pt07_matrix<-as(as.matrix(GetAssayData(Pt07.filter,slot = "counts")), 'sparseMatrix')
Pt07_feature<-data.frame(gene_id=rownames(Pt07_matrix),gene_short_name=rownames(Pt07_matrix))
rownames(Pt07_feature)<-rownames(Pt07_matrix)
Pt07_fd<-new("AnnotatedDataFrame", data = Pt07_feature)
sample_ann<-Pt07.filter@meta.data
rownames(sample_ann)<-colnames(Pt07_matrix)
Pt07_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt07.cds<-newCellDataSet(Pt07_matrix,phenoData =Pt07_pd,featureData =Pt07_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt07.cds <- estimateSizeFactors(Pt07.cds)
Pt07.cds <- estimateDispersions(Pt07.cds)
Pt07.cds<- detectGenes(Pt07.cds, min_expr = 0.1)
Pt07.cds<- setOrderingFilter(Pt07.cds,ordering_genes = EMT)
Pt07.cds <- reduceDimension(
  Pt07.cds,
  max_components = 2, 
  method = 'DDRTree')

Pt07.cds <- orderCells(Pt07.cds)
head(pData(Pt07.cds))
pData(Pt07.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt07.cds <- orderCells(Pt07.cds, root_state = 11)
pData(Pt07.cds)$PseudoEMT<-pData(Pt07.cds)$Pseudotime

plot_cell_trajectory(Pt07.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt07_EMTscore.pdf",width = 8,height = 7,units = "in")

plot_cell_trajectory(Pt07.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt07_PseudoEMT.pdf",width = 10,height = 8,units = "in")


diff_test_res.Pt07.total <- differentialGeneTest(Pt07.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt07.filter<-subset(diff_test_res.Pt07.total, qval < 0.05)
sig_gene_names.Pt07 <- row.names(subset(diff_test_res.Pt07.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt07.cds[sig_gene_names.Pt07,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt07.cds[sig_gene_names.Pt07,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt07.filter=left_join(tep, diff_test_res.Pt07.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt07.filter$Cluster<-factor(diff_test_res.Pt07.filter$Cluster,levels = c(1,5,2,4,6,3), ordered=T) 
b<-diff_test_res.Pt07.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt07.filter$Cluster))), 
         main='Pt07 PseudoEMT DEGs',
         filename="Pt07 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()




msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(1,5,2,4,6,3)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt07=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt07)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt07=diff_test_res.Pt07.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt07$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt07'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt07=rbind(re.Pt07, res@result)
}

Pt07.order=data.frame(Cluster=c(1,5,2,4,6,3), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt07.filter$Cluster=as.numeric(diff_test_res.Pt07.filter$Cluster)
diff_test_res.Pt07.filter=diff_test_res.Pt07.filter %>% left_join(Pt07.order,by='Cluster')
diff_test_res.Pt07.filter$Pt='Pt07'




####Pt08######
Pt08<-Tumor %>% subset(patient%in%'Pt08')
Pt08_matrix<-as(as.matrix(GetAssayData(Pt08,slot = "counts")), 'sparseMatrix')
Pt08_feature<-data.frame(gene_id=rownames(Pt08_matrix),gene_short_name=rownames(Pt08_matrix))
rownames(Pt08_feature)<-rownames(Pt08_matrix)
Pt08_fd<-new("AnnotatedDataFrame", data = Pt08_feature)
sample_ann<-Pt08@meta.data
rownames(sample_ann)<-colnames(Pt08_matrix)
Pt08_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt08.cds<-newCellDataSet(Pt08_matrix,phenoData =Pt08_pd,featureData =Pt08_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt08.cds <- estimateSizeFactors(Pt08.cds)
Pt08.cds <- estimateDispersions(Pt08.cds)
Pt08.cds<- detectGenes(Pt08.cds, min_expr = 0.1)
Pt08.cds<- setOrderingFilter(Pt08.cds,ordering_genes = EMT)
Pt08.cds <- reduceDimension(
  Pt08.cds,
  max_components = 2, 
  method = 'DDRTree')
Pt08.cds <- orderCells(Pt08.cds)
head(pData(Pt08.cds))
pData(Pt08.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt08.cds <- orderCells(Pt08.cds, root_state = 4)
pData(Pt08.cds)$PseudoEMT<-pData(Pt08.cds)$Pseudotime

plot_cell_trajectory(Pt08.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt08_EMTscore.pdf",width = 8,height = 7,units = "in")

plot_cell_trajectory(Pt08.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt08_PseudoEMT.pdf",width = 10,height = 8,units = "in")


diff_test_res.Pt08.total <- differentialGeneTest(Pt08.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt08.filter<-subset(diff_test_res.Pt08.total, qval < 0.05)
sig_gene_names.Pt08 <- row.names(subset(diff_test_res.Pt08.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt08.cds[sig_gene_names.Pt08,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt08.cds[sig_gene_names.Pt08,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt08.filter=left_join(tep, diff_test_res.Pt08.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt08.filter$Cluster<-factor(diff_test_res.Pt08.filter$Cluster,levels = c(3,4,1,6,2,5), ordered=T) 
b<-diff_test_res.Pt08.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt08.filter$Cluster))), 
         main='Pt08 PseudoEMT DEGs',
         filename="Pt08 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()





msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(3,4,1,6,2,5)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt08=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt08)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt08=diff_test_res.Pt08.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt08$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt08'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt08=rbind(re.Pt08, res@result)
}

Pt08.order=data.frame(Cluster=c(3,4,1,6,2,5), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt08.filter$Cluster=as.numeric(diff_test_res.Pt08.filter$Cluster)
diff_test_res.Pt08.filter=diff_test_res.Pt08.filter %>% left_join(Pt08.order,by='Cluster')
diff_test_res.Pt08.filter$Pt='Pt08'



###merge Pt05-08#####
Hallmark.part2=re.Pt05 %>% rbind(re.Pt06) %>% rbind(re.Pt07) %>% rbind(re.Pt08)
#saveRDS(Hallmark.part2,"Hallmark_part2_1127.rds")

diff_test_res.part2=diff_test_res.Pt05.filter %>% 
  rbind(diff_test_res.Pt06.filter) %>% 
  rbind(diff_test_res.Pt07.filter) %>% 
  rbind(diff_test_res.Pt08.filter) 
#saveRDS(diff_test_res.part2,"diff_test_res_part2_1127.rds")

pData.part2=pData(Pt05.cds) %>% rbind(pData(Pt06.cds)) %>% rbind(pData(Pt07.cds)) %>% rbind(pData(Pt08.cds))
#saveRDS(pData.part2, "pData_part2_1129.rds")






####Pt10######
Pt10<-Tumor %>% subset(patient%in%'Pt10')
Pt10_matrix<-as(as.matrix(GetAssayData(Pt10,slot = "counts")), 'sparseMatrix')
Pt10_feature<-data.frame(gene_id=rownames(Pt10_matrix),gene_short_name=rownames(Pt10_matrix))
rownames(Pt10_feature)<-rownames(Pt10_matrix)
Pt10_fd<-new("AnnotatedDataFrame", data = Pt10_feature)
sample_ann<-Pt10@meta.data
rownames(sample_ann)<-colnames(Pt10_matrix)
Pt10_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt10.cds<-newCellDataSet(Pt10_matrix,phenoData =Pt10_pd,featureData =Pt10_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt10.cds <- estimateSizeFactors(Pt10.cds)
Pt10.cds <- estimateDispersions(Pt10.cds)
Pt10.cds<- detectGenes(Pt10.cds, min_expr = 0.1)
Pt10.cds<- setOrderingFilter(Pt10.cds,ordering_genes = EMT)
Pt10.cds <- reduceDimension(
  Pt10.cds,
  max_components = 2,
  method = 'DDRTree')
Pt10.cds <- orderCells(Pt10.cds)
plot_cell_trajectory(Pt10.cds,cell_size = 0.5,color_by = "Pseudotime")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(3,'mm'))+labs(title = 'EMT')+coord_fixed(ratio = 2)+NoLegend() 

dim(Pt10.cds@reducedDimS)
plot(Pt10.cds@reducedDimS[1,],Pt10.cds@reducedDimS[2,])
plot(Pt10.cds@reducedDimW[1,],Pt10.cds@reducedDimW[2,])
sort(Pt10.cds@reducedDimS[2,]) %>% tail(10)
names(sort(Pt10.cds@reducedDimS[2,]) %>% tail(10))

Pt10.filter<-subset(Pt10,cells=names(sort(Pt10.cds@reducedDimS[2,]) %>% tail(5)),invert=T)

Pt10_matrix<-as(as.matrix(GetAssayData(Pt10.filter,slot = "counts")), 'sparseMatrix')
Pt10_feature<-data.frame(gene_id=rownames(Pt10_matrix),gene_short_name=rownames(Pt10_matrix))
rownames(Pt10_feature)<-rownames(Pt10_matrix)
Pt10_fd<-new("AnnotatedDataFrame", data = Pt10_feature)
sample_ann<-Pt10.filter@meta.data
rownames(sample_ann)<-colnames(Pt10_matrix)
Pt10_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt10.cds<-newCellDataSet(Pt10_matrix,phenoData =Pt10_pd,featureData =Pt10_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt10.cds <- estimateSizeFactors(Pt10.cds)
Pt10.cds <- estimateDispersions(Pt10.cds)
Pt10.cds<- detectGenes(Pt10.cds, min_expr = 0.1)
Pt10.cds<- setOrderingFilter(Pt10.cds,ordering_genes = EMT)
Pt10.cds <- reduceDimension(
  Pt10.cds,
  max_components = 2, 
  method = 'DDRTree')

Pt10.cds <- orderCells(Pt10.cds)

pData(Pt10.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt10.cds <- orderCells(Pt10.cds, root_state = 4)
pData(Pt10.cds)$PseudoEMT<-pData(Pt10.cds)$Pseudotime

plot_cell_trajectory(Pt10.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt10_EMTscore.pdf",width = 8,height = 7,units = "in")

plot_cell_trajectory(Pt10.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt10_PseudoEMT.pdf",width = 10,height = 8,units = "in")


diff_test_res.Pt10.total <- differentialGeneTest(Pt10.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt10.filter<-subset(diff_test_res.Pt10.total, qval < 0.05)
sig_gene_names.Pt10 <- row.names(subset(diff_test_res.Pt10.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt10.cds[sig_gene_names.Pt10,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt10.cds[sig_gene_names.Pt10,],  
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt10.filter=left_join(tep, diff_test_res.Pt10.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt10.filter$Cluster<-factor(diff_test_res.Pt10.filter$Cluster,levels = c(1,2,4,6,3,5), ordered=T) 
b<-diff_test_res.Pt10.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt10.filter$Cluster))), 
         main='Pt10 PseudoEMT DEGs',
         filename="Pt10 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()


msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(1,2,4,6,3,5)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt10=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt10)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt10=diff_test_res.Pt10.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt10$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt10'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt10=rbind(re.Pt10, res@result)
}

Pt10.order=data.frame(Cluster=c(1,2,4,6,3,5), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt10.filter$Cluster=as.numeric(diff_test_res.Pt10.filter$Cluster)
diff_test_res.Pt10.filter=diff_test_res.Pt10.filter %>% left_join(Pt10.order,by='Cluster')
diff_test_res.Pt10.filter$Pt='Pt10'




####Pt11######
Pt11<-Tumor %>% subset(patient%in%'Pt11')
Pt11_matrix<-as(as.matrix(GetAssayData(Pt11,slot = "counts")), 'sparseMatrix')
Pt11_feature<-data.frame(gene_id=rownames(Pt11_matrix),gene_short_name=rownames(Pt11_matrix))
rownames(Pt11_feature)<-rownames(Pt11_matrix)
Pt11_fd<-new("AnnotatedDataFrame", data = Pt11_feature)
sample_ann<-Pt11@meta.data
rownames(sample_ann)<-colnames(Pt11_matrix)
Pt11_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt11.cds<-newCellDataSet(Pt11_matrix,phenoData =Pt11_pd,featureData =Pt11_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt11.cds <- estimateSizeFactors(Pt11.cds)
Pt11.cds <- estimateDispersions(Pt11.cds)
Pt11.cds<- detectGenes(Pt11.cds, min_expr = 0.1)
Pt11.cds<- setOrderingFilter(Pt11.cds,ordering_genes = EMT)
Pt11.cds <- reduceDimension(
  Pt11.cds,
  max_components = 2,
  method = 'DDRTree')
Pt11.cds <- orderCells(Pt11.cds)
plot_cell_trajectory(Pt11.cds,cell_size = 0.5,color_by = "Pseudotime")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(3,'mm'))+labs(title = 'EMT')+coord_fixed(ratio = 2)+NoLegend() #+ylim(c(-1,2))

dim(Pt11.cds@reducedDimS)
plot(Pt11.cds@reducedDimS[1,],Pt11.cds@reducedDimS[2,])
plot(Pt11.cds@reducedDimW[1,],Pt11.cds@reducedDimW[2,])
sort(Pt11.cds@reducedDimS[2,]) %>% tail(10)
names(sort(Pt11.cds@reducedDimS[2,]) %>% tail(10))

Pt11.filter<-subset(Pt11,cells=names(sort(Pt11.cds@reducedDimS[2,]) %>% tail(2)),invert=T)

Pt11_matrix<-as(as.matrix(GetAssayData(Pt11.filter,slot = "counts")), 'sparseMatrix')
Pt11_feature<-data.frame(gene_id=rownames(Pt11_matrix),gene_short_name=rownames(Pt11_matrix))
rownames(Pt11_feature)<-rownames(Pt11_matrix)
Pt11_fd<-new("AnnotatedDataFrame", data = Pt11_feature)
sample_ann<-Pt11.filter@meta.data
rownames(sample_ann)<-colnames(Pt11_matrix)
Pt11_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt11.cds<-newCellDataSet(Pt11_matrix,phenoData =Pt11_pd,featureData =Pt11_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt11.cds <- estimateSizeFactors(Pt11.cds)
Pt11.cds <- estimateDispersions(Pt11.cds)
Pt11.cds<- detectGenes(Pt11.cds, min_expr = 0.1)
Pt11.cds<- setOrderingFilter(Pt11.cds,ordering_genes = EMT)
Pt11.cds <- reduceDimension(
  Pt11.cds,
  max_components = 2, 
  method = 'DDRTree')

Pt11.cds <- orderCells(Pt11.cds)

pData(Pt11.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt11.cds <- orderCells(Pt11.cds, root_state = 6)
pData(Pt11.cds)$PseudoEMT<-pData(Pt11.cds)$Pseudotime

plot_cell_trajectory(Pt11.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt11_EMTscore.pdf",width = 8,height = 7,units = "in")


plot_cell_trajectory(Pt11.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt11_PseudoEMT.pdf",width = 10,height = 8,units = "in")



diff_test_res.Pt11.total <- differentialGeneTest(Pt11.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt11.filter<-subset(diff_test_res.Pt11.total, qval < 0.05)
sig_gene_names.Pt11 <- row.names(subset(diff_test_res.Pt11.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt11.cds[sig_gene_names.Pt11,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt11.cds[sig_gene_names.Pt11,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt11.filter=left_join(tep, diff_test_res.Pt11.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt11.filter$Cluster<-factor(diff_test_res.Pt11.filter$Cluster,levels = c(3,5,1,2,6,4), ordered=T) 
b<-diff_test_res.Pt11.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt11.filter$Cluster))), 
         main='Pt11 PseudoEMT DEGs',
         filename="Pt11 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()


msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(3,5,1,2,6,4)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt11=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt11)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt11=diff_test_res.Pt11.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt11$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt11'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt11=rbind(re.Pt11, res@result)
}

Pt11.order=data.frame(Cluster=c(3,5,1,2,6,4), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt11.filter$Cluster=as.numeric(diff_test_res.Pt11.filter$Cluster)
diff_test_res.Pt11.filter=diff_test_res.Pt11.filter %>% left_join(Pt11.order,by='Cluster')
diff_test_res.Pt11.filter$Pt='Pt11'



####Pt15######
Pt15<-Tumor %>% subset(patient%in%'Pt15')
Pt15_matrix<-as(as.matrix(GetAssayData(Pt15,slot = "counts")), 'sparseMatrix')
Pt15_feature<-data.frame(gene_id=rownames(Pt15_matrix),gene_short_name=rownames(Pt15_matrix))
rownames(Pt15_feature)<-rownames(Pt15_matrix)
Pt15_fd<-new("AnnotatedDataFrame", data = Pt15_feature)
sample_ann<-Pt15@meta.data
rownames(sample_ann)<-colnames(Pt15_matrix)
Pt15_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt15.cds<-newCellDataSet(Pt15_matrix,phenoData =Pt15_pd,featureData =Pt15_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt15.cds <- estimateSizeFactors(Pt15.cds)
Pt15.cds <- estimateDispersions(Pt15.cds)
Pt15.cds<- detectGenes(Pt15.cds, min_expr = 0.1)
Pt15.cds<- setOrderingFilter(Pt15.cds,ordering_genes = EMT)
Pt15.cds <- reduceDimension(
  Pt15.cds,
  max_components = 2, 
  method = 'DDRTree')
Pt15.cds <- orderCells(Pt15.cds)

pData(Pt15.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score))
Pt15.cds <- orderCells(Pt15.cds, root_state = 1)
pData(Pt15.cds)$PseudoEMT<-pData(Pt15.cds)$Pseudotime

plot_cell_trajectory(Pt15.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt15_EMTscore.pdf",width = 8,height = 7,units = "in")

plot_cell_trajectory(Pt15.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt15_PseudoEMT.pdf",width = 10,height = 8,units = "in")


diff_test_res.Pt15.total <- differentialGeneTest(Pt15.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt15.filter<-subset(diff_test_res.Pt15.total, qval < 0.05)
sig_gene_names.Pt15 <- row.names(subset(diff_test_res.Pt15.total, qval < 0.05)) 


t <- plot_pseudotime_heatmap(Pt15.cds[sig_gene_names.Pt15,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt15.cds[sig_gene_names.Pt15,],   
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt15.filter=left_join(tep, diff_test_res.Pt15.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt15.filter$Cluster<-factor(diff_test_res.Pt15.filter$Cluster,levels = c(2,5,6,4,1,3), ordered=T) 
b<-diff_test_res.Pt15.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt15.filter$Cluster))), 
         main='Pt15 PseudoEMT DEGs',
         filename="Pt15 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()


msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(2,5,6,4,1,3)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt15=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt15)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt15=diff_test_res.Pt15.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt15$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt15'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt15=rbind(re.Pt15, res@result)
}

Pt15.order=data.frame(Cluster=c(2,5,6,4,1,3), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt15.filter$Cluster=as.numeric(diff_test_res.Pt15.filter$Cluster)
diff_test_res.Pt15.filter=diff_test_res.Pt15.filter %>% left_join(Pt15.order,by='Cluster')
diff_test_res.Pt15.filter$Pt='Pt15'



####Pt17######
Pt17<-Tumor %>% subset(patient%in%'Pt17')
Pt17_matrix<-as(as.matrix(GetAssayData(Pt17,slot = "counts")), 'sparseMatrix')
Pt17_feature<-data.frame(gene_id=rownames(Pt17_matrix),gene_short_name=rownames(Pt17_matrix))
rownames(Pt17_feature)<-rownames(Pt17_matrix)
Pt17_fd<-new("AnnotatedDataFrame", data = Pt17_feature)
sample_ann<-Pt17@meta.data
rownames(sample_ann)<-colnames(Pt17_matrix)
Pt17_pd<-new("AnnotatedDataFrame", data =sample_ann)
Pt17.cds<-newCellDataSet(Pt17_matrix,phenoData =Pt17_pd,featureData =Pt17_fd,lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size())
Pt17.cds <- estimateSizeFactors(Pt17.cds)
Pt17.cds <- estimateDispersions(Pt17.cds)
Pt17.cds<- detectGenes(Pt17.cds, min_expr = 0.1)
Pt17.cds<- setOrderingFilter(Pt17.cds,ordering_genes = EMT)
Pt17.cds <- reduceDimension(
  Pt17.cds,
  max_components = 2, 
  method = 'DDRTree')
Pt17.cds <- orderCells(Pt17.cds)

pData(Pt17.cds) %>% group_by(State) %>% dplyr::summarise(n=mean(EMT_score)) 
Pt17.cds <- orderCells(Pt17.cds, root_state = 6)
pData(Pt17.cds)$PseudoEMT<-pData(Pt17.cds)$Pseudotime

plot_cell_trajectory(Pt17.cds,cell_size = 0.5,color_by = "EMT_score")+coord_fixed(ratio = 1)+
  scale_color_gradientn(colours = c('blue','yellow','red'))+theme(legend.key.width = unit(5,'mm'),legend.key.height = unit(8,'mm'),legend.position = 'right')
ggsave("Pt17_EMTscore.pdf",width = 8,height = 7,units = "in")


plot_cell_trajectory(Pt17.cds,cell_size = 0.5,color_by = "PseudoEMT")+facet_wrap(~lesions)+scale_color_gradientn(colours = c('blue','yellow','red'))+
  theme(legend.position = 'right',legend.key.size = unit(5,'mm'))+coord_fixed(ratio = 1) 
ggsave("Pt17_PseudoEMT.pdf",width = 10,height = 8,units = "in")



diff_test_res.Pt17.total <- differentialGeneTest(Pt17.cds, 
                                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res.Pt17.filter<-subset(diff_test_res.Pt17.total, qval < 0.05)
sig_gene_names.Pt17 <- row.names(subset(diff_test_res.Pt17.total, qval < 0.05))


t <- plot_pseudotime_heatmap(Pt17.cds[sig_gene_names.Pt17,],
                             num_clusters = 6,  
                             cores = 1,
                             show_rownames = T,
                             return_heatmap=T)

a<-pseudotime(Pt17.cds[sig_gene_names.Pt17,],  
              num_clusters = 6,
              cores = 1,
              show_rownames = T,
              return_heatmap=T)

tep <- as.data.frame(cutree(t$tree_row, k=6)) 
colnames(tep) <- "Cluster"
tep$gene_id <- rownames(tep)
diff_test_res.Pt17.filter=left_join(tep, diff_test_res.Pt17.filter,by='gene_id')
table(tep$Cluster)

diff_test_res.Pt17.filter$Cluster<-factor(diff_test_res.Pt17.filter$Cluster,levels = c(5,6,3,1,2,4), ordered=T) 
b<-diff_test_res.Pt17.filter %>% arrange(Cluster)
pheatmap(a[b$gene_id,], 
         color = colorRampPalette(colors = c('blue','white','red'))(100),
         cluster_cols = F, cluster_rows = F, 
         show_rownames = F,show_colnames = F,
         border_color = NA, gaps_row = c(cumsum(table(diff_test_res.Pt17.filter$Cluster))), 
         main='Pt17 PseudoEMT DEGs',
         filename="Pt17 PseudoEMT DEGs.pdf",width=4,height = 7)
dev.new()


msigdb = msigdf.human %>% select(geneset,symbol) %>% as.data.frame() 
EMP_order=c(5,6,3,1,2,4)
pseudoEMT_order=c('E','H1','H2','H3','H4','M')
re.Pt17=as.data.frame(matrix(nrow=0,ncol = 12))
colnames(re.Pt17)=c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count","Filter","Pt","pseudoEMT")

for (i in 1:6) {
  temp.Pt17=diff_test_res.Pt17.filter %>% filter(Cluster==EMP_order[i])
  res=enricher(temp.Pt17$gene_id, TERM2GENE = msigdb) 
  res@result$Filter=word(res@result$ID,1,sep = '_') 
  res@result$Pt='Pt17'
  res@result$pseudoEMT=pseudoEMT_order[i]
  res@result=res@result %>% filter(Filter=='HALLMARK', qvalue < 0.05) 
  res@result=res@result[1:10,]  
  res@result=na.omit(res@result)  
  rownames(res@result)=NULL
  re.Pt17=rbind(re.Pt17, res@result)
}

Pt17.order=data.frame(Cluster=c(5,6,3,1,2,4), pseudoEMT_order=c('E','H1','H2','H3','H4','M'))
diff_test_res.Pt17.filter$Cluster=as.numeric(diff_test_res.Pt17.filter$Cluster)
diff_test_res.Pt17.filter=diff_test_res.Pt17.filter %>% left_join(Pt17.order,by='Cluster')
diff_test_res.Pt17.filter$Pt='Pt17'


##merge Pt10-11-15-17#####
Hallmark.part3=re.Pt10 %>% rbind(re.Pt11) %>% rbind(re.Pt15) %>% rbind(re.Pt17)
#saveRDS(Hallmark.part3,"Hallmark_part3_1127.rds")

diff_test_res.part3=diff_test_res.Pt10.filter %>% 
  rbind(diff_test_res.Pt11.filter) %>% 
  rbind(diff_test_res.Pt15.filter) %>% 
  rbind(diff_test_res.Pt17.filter) 
#saveRDS(diff_test_res.part3,"Tumor_filter/diff_test_res_part3_1127.rds")


pData.part3=pData(Pt10.cds) %>% rbind(pData(Pt11.cds)) %>% rbind(pData(Pt15.cds)) %>% rbind(pData(Pt17.cds))
#saveRDS(pData.part3, "pData_part3_1129.rds")







#Fig 5E####
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggsci)
library(ggsignif)

pData=rbind(pData.part1) %>% rbind(pData.part2) %>% rbind(pData.part3)

pData.copy=pData
pData.new=as.data.frame(matrix(nrow=0, ncol = 23))
colnames(pData.new)=c(colnames(pData.copy),'PseudoEMT.new')
pt.list=unique(pData$patient)
for (i in pt.list) {
  pD.tmp=subset(pData, patient%in%i)
  pD.tmp$PseudoEMT.new=pD.tmp$PseudoEMT/max(pD.tmp$PseudoEMT)*100   
  pData.new=rbind(pData.new, pD.tmp) 
}
pData.new$barcode=rownames(pData.new)


head(pData.new)
dim(pData.new) 

Tumor = readRDS("tumors_filter_by_cnv_man.rds")

IG_name<-rownames(Tumor)[grep("^IG[HLK]",rownames(Tumor))]
HSP_name<-rownames(Tumor)[grep("^HSP",rownames(Tumor))]
keep_name<-setdiff(rownames(Tumor), c(IG_name,"JCHAIN")) %>% setdiff(HSP_name)
Tumor<-subset(Tumor, features=keep_name)

excludeGene<-read.table("excludeGene.txt",header=TRUE)
Tumor<-Tumor[!rownames(Tumor) %in% excludeGene$Gene]


Tumor=Tumor %>% subset(cells=rownames(pData.new))

Tumor$PseudoEMT=pData.new$PseudoEMT.new
Tumor$EMT_score_old=pData.new$EMT_score

metadata=Tumor@meta.data
head(metadata)
dim(metadata)
pt.list=unique(metadata$patient)

EMT<-read.delim('EMT.txt')
EMT<-EMT[-1,]
EMT.list<-list(EMT)

seurat.tumor=Tumor %>% subset(patient%in%'Pt01')
seurat.tumor=AddModuleScore(object = seurat.tumor,features = EMT.list, name='EMT_score')
for (i in pt.list[2:13]) {
  seurat.tmp=Tumor %>% subset(patient%in%i)
  seurat.tmp=AddModuleScore(object = seurat.tmp,features = EMT.list, name='EMT_score')
  seurat.tumor=merge(seurat.tumor,seurat.tmp)
}
seurat.tumor

head(seurat.tumor)
colnames(seurat.tumor@meta.data)[17]='EMT_score_new'

#计算前后的各20%细胞

metadata=seurat.tumor@meta.data

EMT_compare_E=c()
EMT_compare_H1=c()
EMT_compare_H2=c()
EMT_compare_H3=c()
EMT_compare_H4=c()
EMT_compare_M=c()


for (i in pt.list) {
  pD.tmp=subset(metadata, patient%in%i)
  pD.E =pD.tmp%>% arrange(PseudoEMT) %>% subset(PseudoEMT<unname(quantile(pD.tmp$PseudoEMT,1/6))) %>% summarise(EMT_score=mean(EMT_score_new))
  pD.H1 =pD.tmp%>% arrange(PseudoEMT) %>% subset(PseudoEMT<unname(quantile(pD.tmp$PseudoEMT,2/6))) %>% subset(PseudoEMT>unname(quantile(pD.tmp$PseudoEMT,1/6))) %>% summarise(EMT_score=mean(EMT_score_new))
  pD.H2 =pD.tmp%>% arrange(PseudoEMT) %>% subset(PseudoEMT<unname(quantile(pD.tmp$PseudoEMT,3/6))) %>% subset(PseudoEMT>unname(quantile(pD.tmp$PseudoEMT,2/6))) %>% summarise(EMT_score=mean(EMT_score_new))
  pD.H3 =pD.tmp%>% arrange(PseudoEMT) %>% subset(PseudoEMT<unname(quantile(pD.tmp$PseudoEMT,4/6))) %>% subset(PseudoEMT>unname(quantile(pD.tmp$PseudoEMT,3/6))) %>% summarise(EMT_score=mean(EMT_score_new))
  pD.H4 =pD.tmp%>% arrange(PseudoEMT) %>% subset(PseudoEMT<unname(quantile(pD.tmp$PseudoEMT,5/6))) %>% subset(PseudoEMT>unname(quantile(pD.tmp$PseudoEMT,4/6))) %>% summarise(EMT_score=mean(EMT_score_new))
  pD.M =pD.tmp%>% arrange(PseudoEMT) %>% subset(PseudoEMT>unname(quantile(pD.tmp$PseudoEMT,5/6))) %>% summarise(EMT_score=mean(EMT_score_new))
  EMT_compare_E=rbind(EMT_compare_E, pD.E)
  EMT_compare_H1=rbind(EMT_compare_H1, pD.H1)
  EMT_compare_H2=rbind(EMT_compare_H2, pD.H2)
  EMT_compare_H3=rbind(EMT_compare_H3, pD.H3)
  EMT_compare_H4=rbind(EMT_compare_H4, pD.H4)
  EMT_compare_M=rbind(EMT_compare_M, pD.M)
}


EMT_compare=data.frame(
  EMT_score=c(EMT_compare_E$EMT_score,EMT_compare_H1$EMT_score,EMT_compare_H2$EMT_score,EMT_compare_H3$EMT_score,EMT_compare_H4$EMT_score,EMT_compare_M$EMT_score),
  Group=c(rep('E',13),rep('H1',13),rep('H2',13),rep('H3',13),rep('H4',13),rep('M',13)))


EMT_compare$Phase[EMT_compare$Group=='E']=1
EMT_compare$Phase[EMT_compare$Group=='H1']=2
EMT_compare$Phase[EMT_compare$Group=='H2']=3
EMT_compare$Phase[EMT_compare$Group=='H3']=4
EMT_compare$Phase[EMT_compare$Group=='H4']=5
EMT_compare$Phase[EMT_compare$Group=='M']=6


ggplot(EMT_compare, aes(Phase, EMT_score, color=as.factor(Phase)))+
  geom_boxplot(outlier.fill = "white",outlier.color = "white")+
  geom_jitter(position=position_jitter(width = 0.3, height=0),alpha=1,size=1.5)+
  geom_smooth(method = "lm", se=T, color="black", formula = y~x, size=0.75)+
  ggpubr::stat_cor(data=EMT_compare, mapping = aes(Phase,EMT_score), method='pearson', size=3, inherit.aes = F)+
  theme_bw()+
  scale_color_d3("category20")+
  scale_x_continuous(breaks = seq(1,6, by=1))+
  theme(axis.text=element_text(face = "bold"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))+
  labs(title = "EMT score from E to M phase")
ggsave("Fig 5E.pdf",width = 5,height = 3.8,units = "in")




##Fig 5F#####
library(reshape2)
library(ggplot2)
library(ggpubr)
library(stringr)
Hallmark.res=Hallmark.part1 %>% rbind(Hallmark.part2) %>% rbind(Hallmark.part3)
rm(list = c('Hallmark.part1','Hallmark.part2','Hallmark.part3'))


Hallmark.times=Hallmark.res %>% group_by(pseudoEMT, Description) %>%summarise(n=n())
Hallmark.times=Hallmark.times %>% dcast(Description~pseudoEMT)
Hallmark.times[is.na(Hallmark.times)]=0
row.names(Hallmark.times)=Hallmark.times$Description
Hallmark.times=Hallmark.times[,-1]

Hallmark.res.copy=Hallmark.res
for (i in 1:dim(Hallmark.res.copy)[1]) {
  Hallmark.res.copy$GeneRatio.num[i]= as.numeric(strsplit(Hallmark.res.copy$GeneRatio,'/')[[i]][1]) / as.numeric(strsplit(Hallmark.res.copy$GeneRatio,'/')[[i]][2])
}

Hallmark.times=Hallmark.res.copy %>% group_by(pseudoEMT, Description) %>% summarise(sum=sum(GeneRatio.num,na.rm=T)) 

Hallmark.times=Hallmark.times %>% dcast(Description~pseudoEMT)
Hallmark.times[is.na(Hallmark.times)]=0
row.names(Hallmark.times)=Hallmark.times$Description
Hallmark.times=Hallmark.times[,-1]


Hallmark.times=Hallmark.times[c('HALLMARK_G2M_CHECKPOINT',
                                'HALLMARK_SPERMATOGENESIS',
                                'HALLMARK_E2F_TARGETS',
                                'HALLMARK_MITOTIC_SPINDLE',
                                'HALLMARK_ALLOGRAFT_REJECTION',
                                'HALLMARK_FATTY_ACID_METABOLISM',
                                'HALLMARK_XENOBIOTIC_METABOLISM',
                                'HALLMARK_MYC_TARGETS_V1',
                                'HALLMARK_MYC_TARGETS_V2',
                                'HALLMARK_PEROXISOME',
                                'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                                'HALLMARK_UNFOLDED_PROTEIN_RESPONSE',
                                'HALLMARK_GLYCOLYSIS',
                                'HALLMARK_DNA_REPAIR',
                                'HALLMARK_ADIPOGENESIS',
                                'HALLMARK_NOTCH_SIGNALING',
                                'HALLMARK_PROTEIN_SECRETION',
                                'HALLMARK_MYOGENESIS',
                                'HALLMARK_PI3K_AKT_MTOR_SIGNALING',
                                'HALLMARK_P53_PATHWAY',
                                'HALLMARK_MTORC1_SIGNALING',
                                'HALLMARK_WNT_BETA_CATENIN_SIGNALING',
                                'HALLMARK_ANGIOGENESIS',
                                'HALLMARK_KRAS_SIGNALING_UP',
                                'HALLMARK_TGF_BETA_SIGNALING',
                                'HALLMARK_IL6_JAK_STAT3_SIGNALING',
                                'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY',
                                'HALLMARK_HYPOXIA',
                                'HALLMARK_APICAL_JUNCTION',
                                'HALLMARK_IL2_STAT5_SIGNALING',
                                'HALLMARK_APOPTOSIS',
                                'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
                                'HALLMARK_INTERFERON_ALPHA_RESPONSE',
                                'HALLMARK_INTERFERON_GAMMA_RESPONSE',
                                'HALLMARK_COMPLEMENT',
                                'HALLMARK_INFLAMMATORY_RESPONSE',
                                'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION'
),]


pheatmap(Hallmark.times, 
         cellwidth = 20,cellheight = 8,
         cluster_cols = F, cluster_rows = F,
         border_color = 'white',
         scale = 'row', 
         fontsize = 9,
         treeheight_row = 20,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         angle_col = 0,
         cutree_rows = 6,
         main='',
         gaps_row = c(5,10,17,22,27))


rownames(Hallmark.times)=rownames(Hallmark.times) %>% str_sub(10) 

pheatmap(Hallmark.times, 
         cellwidth = 20,cellheight = 8,
         cluster_cols = F, cluster_rows = F,
         border_color = 'white',
         scale = 'row', 
         fontsize = 7,
         treeheight_row = 20,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         angle_col = 0,
         cutree_rows = 6,
         main='',
         gaps_row = c(5,10,17,22,27),
         filename="Fig 5F.pdf",width=5,height = 5)
dev.new()



#Fig S5F######
library(ggridges)
library(ggplot2)
library(ggsci)
library(dplyr)
pData=rbind(pData.part1) %>% rbind(pData.part2) %>% rbind(pData.part3)
rm(list = c('pData.part1','pData.part2','pData.part3'))
head(pData)


pData.new=as.data.frame(matrix(nrow=0, ncol = 23))
colnames(pData.new)=c(colnames(pData),'PseudoEMT.new')
pt.list=unique(pData$patient)
for (i in pt.list) {
  pD.tmp=subset(pData, patient%in%i)
  pD.tmp$PseudoEMT.new=pD.tmp$PseudoEMT/max(pD.tmp$PseudoEMT)*100 
  pData.new=rbind(pData.new, pD.tmp) 
}
pData.new$barcode=rownames(pData.new)


pData.new$lesions=factor(pData.new$lesions,levels = c('Liver_meta','Colon_tumor'),ordered = T)
ggplot(pData.new,aes(x=PseudoEMT.new, y=lesions))+
  geom_density_ridges_gradient(aes(fill=..x..))+
  scale_fill_gradientn(colors = c('blue','white','red'),breaks=c(0,50,100))+
  labs(title = '')+theme_bw()+theme(plot.title = element_text(hjust = 0.5))+
  ylab('')+xlab('PseudoEMT')+
  scale_x_continuous(breaks = seq(0,100, by=100/6))+
  guides(fill=guide_colorbar(title='PseudoEMT'))+
  theme(axis.text=element_text(face = "bold",colour = "black"),
        axis.text.x = element_text(angle = 0,vjust=0,hjust=0.5),
        plot.title = element_text(hjust = 0.5,size = 12,face = "bold",colour = "black"))+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
ggsave("Fig S5F.pdf",width = 6,height = 4,units = "in")




#Fig S5G######

pData.new$phase[pData.new$PseudoEMT.new<=100/6]='E'
pData.new$phase[pData.new$PseudoEMT.new>100/6 & pData.new$PseudoEMT.new<=200/6]='H1'
pData.new$phase[pData.new$PseudoEMT.new>200/6 & pData.new$PseudoEMT.new<=300/6]='H2'
pData.new$phase[pData.new$PseudoEMT.new>300/6 & pData.new$PseudoEMT.new<=400/6]='H3'
pData.new$phase[pData.new$PseudoEMT.new>400/6 & pData.new$PseudoEMT.new<=500/6]='H4'
pData.new$phase[pData.new$PseudoEMT.new>500/6 & pData.new$PseudoEMT.new<=100]='M'

table(pData.new$lesions,pData.new$phase,pData.new$patient)
table(pData.new$lesions,pData.new$patient)

ratio_data=as.data.frame(table(pData.new$lesions,pData.new$phase,pData.new$patient))
colnames(ratio_data)=c('lesions','phase','patient','Freq')
ratio_data$Pt_lesion=paste(ratio_data$patient,'_',ratio_data$lesions)

ratio_data$Pt_lesion=factor(ratio_data$Pt_lesion,levels=unique(ratio_data$Pt_lesion),ordered = T)
ratio_data$patient=factor(ratio_data$patient,levels=rev(unique(ratio_data$patient)),ordered = T)

ggplot(ratio_data, aes(fill=phase, y=Pt_lesion, x=Freq)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_d3('category20')+theme_light()+
  theme(axis.title =element_text(size = 8),axis.text =element_text(size = 8, color = 'black'))+
  theme(axis.text.x = element_blank(), strip.background = element_blank(),strip.text = element_blank(),
        axis.ticks.x = element_blank())+scale_x_reverse() +facet_wrap(~patient, scales='free_y',ncol=1)+
  ylab('')+xlab('')+labs(title = 'Proportion')+theme(plot.title = element_text(hjust = 0.5))

ggsave("Fig S5G.pdf",width = 6,height = 8,units = "in")





