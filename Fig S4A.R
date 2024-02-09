library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Seurat)
library(ggsci)
library(reshape2)
library(ggpubr)

list_to_vector<-function(ls){
  tmp = c()
  for(i in 1:length(ls)){
    tmp = c(tmp, as.vector(ls[[i]]))
  }
  return(tmp)
}

col_map_func <- function(dat){
  res = list()
  for(n in colnames(dat)){
    tmp_name = sort(unique(dat[, n]))
    tmp_col = pal_igv()(length(tmp_name))
    names(tmp_col) = tmp_name
    res[[n]] = tmp_col
  }
  return(res)
}


infercnv_obj = readRDS('run.final.infercnv_obj')
expr <- infercnv_obj@expr.data
normal_loc <- list_to_vector(infercnv_obj@reference_grouped_cell_indices)
test_loc <- list_to_vector(infercnv_obj@observation_grouped_cell_indices)

anno.df=data.frame(
  CB=c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class=c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)

gn <- rownames(expr)
genefile_dir = "hg38_filtered_order_gene.txt"
geneFile <- read.table(genefile_dir)

rownames(geneFile)=geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr=expr[intersect(gn,geneFile$V1),]
head(sub_geneFile,4)
expr[1:4,1:4]
malignant_tumor = readRDS("tumors_filter_by_cnv_man.rds")
malignant_tumor = subset(malignant_tumor, lesions%in%c('Colon_tumor', 'Liver_meta'))


D1 = expr[,intersect(colnames(malignant_tumor), colnames(expr))]
D1_meta = malignant_tumor@meta.data[colnames(D1), ]
left_anno = D1_meta[, c('patient', 'lesions')]

left_anno_col = col_map_func(left_anno)
clone_order =rownames(left_anno %>% arrange(desc(lesions)))

left_anno = left_anno[clone_order,]

left_anno = rowAnnotation(df = left_anno,
                          col = left_anno_col)
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
bk<-read.table("infercnv.heatmap_thresholds.txt", header=F)[,1]
cols<-colorRampPalette(colors = c("darkblue", "white", "darkred"))(length(bk))
ht = Heatmap(t(D1)[clone_order, ], 
             col=circlize::colorRamp2(bk, cols),
             cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
             column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
             column_gap = unit(2, "mm"),
             heatmap_legend_param = list(
               title = "CNV",
               direction = "vertical",
               title_position = "leftcenter-rot",
               legend_height = unit(3, "cm")),
             top_annotation = top_anno,
             left_annotation = left_anno, 
             row_title = NULL,column_title = NULL)
pdf("Fig S4A.pdf",width = 15,height = 10)
draw(ht, heatmap_legend_side = "right")
dev.off()





