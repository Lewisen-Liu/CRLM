library(ggplot2)
library(dplyr)
library(Seurat)
library(ggsci)
setwd('/data share/data/')
vlnplot_modify <- function(obj, marker_genes, group_by){
  ######
  # @author        : xinwang 
  # @version       : 1.0.0
  # @desc          : vlnplot moodify
  # @param obj     : Seurat object. 
  # @param marker_genes     : a vector of genes.
  # @param group_by     : column name in seurat, such as 'seurat_cluster'.
  # @return        : ggplot object. 
  ######
  gene_exp = FetchData(obj, vars=marker_genes)
  gene_exp$group = obj@meta.data[rownames(gene_exp), group_by]
  
  gene_exp_mean = gene_exp %>%
    group_by(group) %>%
    summarise_all(list(mean)) %>% as.data.frame()
  gene_exp_mean = melt(gene_exp_mean)
  colnames(gene_exp_mean) = c('group1', 'gene', 'mean_exp')
  plot_data = melt(gene_exp)
  colnames(plot_data) = c('group2', 'gene', 'all_exp')
  
  gene_exp_mean$merge_id = paste0(gene_exp_mean$group1, gene_exp_mean$gene)
  plot_data$merge_id = paste0(plot_data$group2, plot_data$gene)
  plot_data = merge(plot_data, gene_exp_mean, by='merge_id')
  
  res=ggplot(plot_data, aes(x=group1, y=all_exp, fill=mean_exp)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE)+
    scale_fill_gradient(low = "#F2E80D", high = "#E65224") +
    scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(rows = vars(gene.x), scales = "free", switch = "y") +
    theme_cowplot(font_size = 12) +
    theme(panel.spacing = unit(0, "lines"),
          panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.y.left = element_text(angle = 0),
          axis.text.x = element_text(angle=60, vjust=1,hjust=1))
  return(res)
}
library(Seurat)
library(iTALK)
library(monocle)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggtext)
vlnplot_modify_group <- function(obj, marker_genes, group_by1, group_by2){
  ######
  # @author        : xinwang 
  # @version       : 1.0.0
  # @desc          : vlnplot moodify
  # @param obj     : Seurat object. 
  # @param marker_genes     : a vector of genes.
  # @param group_by1     : column name in seurat, such as 'seurat_cluster'.
  # @param group_by2     : column name in seurat, such as 'seurat_cluster'.
  # @return        : ggplot object. 
  ######
  gene_exp = FetchData(obj, vars=marker_genes)
  gene_exp$group = obj@meta.data[rownames(gene_exp), group_by1]
  gene_exp$group2 = obj@meta.data[rownames(gene_exp), group_by2]
  
  gene_exp_mean = gene_exp %>%
    group_by(group, group2) %>%
    summarise_all(list(mean)) %>% as.data.frame()
  
  gene_exp_mean = melt(gene_exp_mean)
  colnames(gene_exp_mean) = c('group_by1', group_by2, 'gene', 'mean_exp')
  plot_data = melt(gene_exp)
  colnames(plot_data) = c('group_by1', group_by2, 'gene', 'all_exp')
  
  plot_data = merge(plot_data, gene_exp_mean)
  
  res=ggplot(plot_data, aes_string(x='group_by1', y='all_exp', fill='mean_exp', color=group_by2)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE)+
    scale_fill_gradient(low = "#F2E80D", high = "#E65224") +
    ggsci::scale_color_igv() +
    scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    facet_grid(rows = vars(gene),cols = vars(group_by1), scales = "free", switch = "y") +
    theme_cowplot(font_size = 12) +
    theme(panel.spacing.x  = unit(0.1, "lines"),
          panel.spacing.y  = unit(0, "lines"),
          #axis.line.x.bottom = element_line(color = 'black'),
          #axis.line.x.top = element_line(color = 'black'),
          #axis.line.y.right = element_line(color = 'black'),
          panel.background = element_rect(fill = NA, color = "black"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.text.y.left = element_text(angle = 0),
          strip.text.x = element_blank(),
          # strip.placement =  = unit(1, "mm"),
          axis.text.x = element_text(angle=60, vjust=1,hjust=1))
  return(res)
}

library(circlize)
wx_LRPlot_extend<-function(data,datatype,genes=NULL,gene_col=NULL,transparency=0.5,link.arr.lwd=1,link.arr.lty=NULL,link.arr.col=NULL,link.arr.width=NULL,
                           link.arr.type=NULL,facing='clockwise',cell_col=NULL,print.cell=TRUE,track.height_1=uh(2,'mm'),track.height_2=uh(12,'mm'),
                           annotation.height_1=0.01,annotation.height_2=0.01,text.vjust = '0.4cm',...){
  cell_group<-unique(c(data$cell_from,data$cell_to))
  genes<-c(structure(data$ligand,names=data$cell_from),structure(data$receptor,names=data$cell_to))
  genes<-genes[!duplicated(paste(names(genes),genes))]
  genes<-genes[order(names(genes))]
  if(is.null(link.arr.lty)){
    if(datatype=='mean count'){
      link.arr.lty='solid'
    }else if(datatype=='DEG'){
      link.arr.lty=structure(ifelse(data$cell_from_logFC==0.0001,'dashed','solid'),names=paste(data$cell_from,data$receptor))
    }else{
      print('invalid datatype')
    }
  }
  if(is.null(link.arr.col)){
    if(datatype=='mean count'){
      data<-data %>% mutate(link_col='black')
    }else if(datatype=='DEG'){
      data<-data %>% mutate(link_col=ifelse(cell_from_logFC==0.0001,ifelse(cell_to_logFC>0,'#d73027','#00ccff'),
                                            ifelse(cell_to_logFC==0.0001,ifelse(cell_from_logFC>0,'#d73027','#00ccff'),
                                                   ifelse(cell_from_logFC>0,ifelse(cell_to_logFC>0,'#d73027','#dfc27d'),
                                                          ifelse(cell_to_logFC>0,'#9933ff','#00ccff')))))
    }else{
      print('invalid datatype')
    }
  }else{
    data$link_col=link.arr.col
  }
  if(is.null(link.arr.type)){
    if(datatype=='mean count'){
      link.arr.type='triangle'
    }else if(datatype=='DEG'){
      link.arr.type=structure(ifelse(data$cell_to_logFC==0.0001,'ellipse','triangle'),names=paste(data$cell_from,data$receptor))
    }else{
      print('invalid datatype')
    }
  }
  if(is.null(gene_col)){
    comm_col<-structure(c('#99ff99','#99ccff','#ff9999','#ffcc99'),names=c('other','cytokine','checkpoint','growth factor'))
    gene_col<-structure(c(comm_col[data$comm_type],rep('#073c53',length(data$receptor))),names=c(data$ligand,data$receptor))
  }
  if(is.null(cell_col)){
    cell_col<-structure(randomColor(count=length(unique(names(genes))),luminosity='dark'),names=unique(names(genes)))
  }
  if(is.null(link.arr.lwd)){
    data<-data %>% mutate(arr_width=1)
  }else if(max(abs(link.arr.lwd))-min(abs(link.arr.lwd))==0 && all(link.arr.lwd!=0.0001)){
    data<-data %>% mutate(arr_width=ifelse(abs(link.arr.lwd<5),abs(link.arr.lwd),5))
  }else{
    data<-data %>% mutate(arr_width=ifelse(link.arr.lwd==0.0001,2,1+5/(max(abs(link.arr.lwd))-min(abs(link.arr.lwd)))*(abs(link.arr.lwd)-min(abs(link.arr.lwd)))))
  }
  if(length(cell_group)!=1){
    gap.degree <- do.call("c", lapply(table(names(genes)), function(i) c(rep(1, i-1), 8)))
  }else{
    gap.degree <- do.call("c", lapply(table(names(genes)), function(i) c(rep(1, i))))
  }
  circos.par(gap.degree = gap.degree)
  if(length(gene_col)==1){
    grid.col=gene_col
  }else{
    grid.col=gene_col[genes]
    names(grid.col)<-paste(names(genes),genes)
  }
  if(is.null(link.arr.width)){
    data<-data %>% mutate(link.arr.width=data$arr_width/10)
  }else if(max(abs(link.arr.width))-min(abs(link.arr.width))==0 && all(link.arr.width!=0.0001)){
    data<-data %>% mutate(link.arr.width=ifelse(abs(link.arr.width)<0.5,abs(link.arr.width),0.5))
  }else{
    data<-data %>% mutate(link.arr.width=ifelse(link.arr.width==0.0001,0.2,(1+5/(max(abs(link.arr.width))-min(abs(link.arr.width)))*(abs(link.arr.width)-min(abs(link.arr.width))))/10))
  }
  ##### by wangxin ####
  df = data.frame(from = paste(data$cell_from,data$ligand),
                  to  = paste(data$cell_to,data$receptor),
                  value = !(data$cell_from_logFC==0 & data$cell_to_logFC==0))
  #####################
  chordDiagram(df, 
               order=paste(names(genes),genes),
               grid.col=grid.col,
               link.visible = df$value, # added by wangxin
               transparency=transparency,
               directional=1,
               direction.type='arrows',
               link.arr.lwd=data$arr_width,
               link.arr.lty=link.arr.lty,
               link.arr.type=link.arr.type,
               link.arr.width=data$link.arr.width,
               link.arr.col=data$link_col,
               col='#00000000',
               annotationTrack=c('grid'),
               preAllocateTracks = list(
                 list(track.height = track.height_1),
                 list(track.height = track.height_2)),
               annotationTrackHeight = c(annotation.height_1,annotation.height_2),...)
  
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = genes[get.cell.meta.data("sector.numeric.index")]
    circos.text(mean(xlim),mean(ylim),sector.index, col = "black", cex = 0.7, facing = facing, niceFacing = TRUE)
  }, bg.border = 0)
  
  if(print.cell){
    for(c in unique(names(genes))) {
      gene = as.character(genes[names(genes) == c])
      highlight.sector(sector.index = paste(c,gene), track.index = 1, col = ifelse(length(cell_col)==1,cell_col,cell_col[c]), text = c, text.vjust = text.vjust, niceFacing = TRUE,lwd=1)
    }
  }
  circos.clear()
}

chunk_vector <- function(x,n){split(x, factor(sort(rank(x)%%n)))}
large_dgCMatrix_load <- function(obj, slot='counts',split_num=3){
  dgCMatrix_mat = GetAssayData(obj, slot = slot)
  if(split_num == 1){
    return(as.matrix(dgCMatrix_mat))
  }
  inter_seq = chunk_vector(1:ncol(dgCMatrix_mat), split_num)
  res = c()
  for(start in 1:length(inter_seq)){
    if(start==1){
      res = as.matrix(dgCMatrix_mat[, inter_seq[[start]]])
    }else{
      res = cbind(res, as.matrix(dgCMatrix_mat[, inter_seq[[start]]]))
    }
  }
  return(as.matrix(res))
}
run_iTALK<-function(obj, groupby){
  data = large_dgCMatrix_load(obj, slot='data',split_num=6)
  data = as.data.frame(t(data))
  data$cell_type = obj@meta.data[, groupby]
  
  highly_exprs_genes<-rawParse(data,top_genes=50,stats='mean')
  comm_list<-c('growth factor','other','cytokine','checkpoint')
  res<-NULL
  for(comm_type in comm_list){
    res_cat<-FindLR(highly_exprs_genes,datatype='mean count',comm_type=comm_type)
    res_cat<-res_cat[order(res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs,decreasing=T),]
    res<-rbind(res,res_cat)
  }
  res<-res[order(res$cell_from_mean_exprs*res$cell_to_mean_exprs,decreasing=T),]
  return(res)
}
run_iTALK_DEG <- function(data1, data2){
  res<-NULL
  comm_list<-c('growth factor','other','cytokine','checkpoint')
  for(comm_type in comm_list){
    res_cat<-FindLR(data1, data2, datatype='DEG', comm_type=comm_type)
    res_cat<-res_cat[order(res_cat$cell_from_logFC*res_cat$cell_to_logFC,decreasing=T),]
    res<-rbind(res,res_cat)
  }
  res<-res[order(res$cell_from_logFC*res$cell_to_logFC,decreasing=T),]
  return(res)
}


num = 30
# ----- signif sub_cluster and tumor cell interaction ---
tumor = subset(readRDS('./Total_20220309.rds'), cluster_name=='Tumor cell')
tumor[['new_idents']] = 'CRC_tumor'

Regev_tumor = readRDS('Regev_CRC_intratumor_tumor.rds')
Regev_tumor[['new_idents']] = 'Regev_tumor'

all_srt = readRDS('Regev_CRC_project.rds')
all_srt[['new_name']] = NA
all_srt@meta.data[rownames(all_srt@meta.data[all_srt$From=='CRC', ]), 'new_name'] = all_srt@meta.data[all_srt$From=='CRC', 'cluster_name']
all_srt@meta.data[rownames(all_srt@meta.data[all_srt$From=='Regev', ]), 'new_name'] = all_srt@meta.data[all_srt$From=='Regev', 'predicted.cluster_name']
Idents(all_srt) = all_srt$new_name

my_th17 = subset(all_srt, new_name=='CD4_IL17A_Th17'&From=='CRC')
my_th17[['new_idents']] = 'CRC_Th17'

regev_th17 = subset(all_srt, new_name=='CD4_IL17A_Th17'&From=='Regev')
regev_th17[['new_idents']] = 'Regev_Th17'


My_th17_tumor = merge(tumor, my_th17)
Idents(My_th17_tumor) = My_th17_tumor$new_idents

Regev_th17_tumor = merge(Regev_tumor, regev_th17)
Idents(Regev_th17_tumor) = Regev_th17_tumor$new_idents

sub_srt =  merge(My_th17_tumor, Regev_th17_tumor)
Idents(sub_srt) = sub_srt$new_idents
process_srt <-function(obj, n.components=2,resolution=0.5){
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj, features = VariableFeatures(object = obj))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 50, verbose = TRUE)
  obj <- FindNeighbors(obj, reduction = "pca",dims = 1:20)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj,dims = 1:20,reduction = "pca", n.components=n.components)
  return(obj)
}
sub_srt = process_srt(sub_srt)

DimPlot(sub_srt, group.by='new_idents')
VlnPlot(sub_srt, features=c('TNFRSF12A', 'TNFSF12'), pt.size=0, group.by = 'new_idents')
receptor_exp = AverageExpression(sub_srt, features = 'TNFRSF12A', group.by = 'orig.ident')$RNA %>% log1p()

Idents(sub_srt) = sub_srt$new_idents
th17_deg = FindMarkers(subset(sub_srt, new_idents%in%c('Regev_Th17', 'CRC_Th17')), 
                       ident.1 = 'CRC_Th17', ident.2 = 'Regev_Th17')
tumor_deg = FindMarkers(subset(sub_srt, new_idents%in%c('Regev_tumor', 'CRC_tumor')), 
                        ident.1 = 'CRC_tumor', ident.2 = 'Regev_tumor')

th17_deg = th17_deg %>% filter(p_val_adj<0.01)
th17_deg$gene = rownames(th17_deg)
th17_deg = th17_deg[, c('gene', 'avg_log2FC','p_val', 'p_val_adj')]
th17_deg$cell_type = 'Th17'
colnames(th17_deg) = c('gene', 'logFC', 'p.value', 'q.value', 'cell_type')

tumor_deg = tumor_deg %>% filter(p_val_adj<0.01)
tumor_deg$gene = rownames(tumor_deg)
tumor_deg = tumor_deg[, c('gene', 'avg_log2FC','p_val', 'p_val_adj')]
tumor_deg$cell_type = 'Tumor'
colnames(tumor_deg) = c('gene', 'logFC', 'p.value', 'q.value', 'cell_type')

saveRDS(tumor_deg, './italk_tumor_deg_cell.rds')
saveRDS(th17_deg, './italk_th17_deg_cell.rds')


tumor_deg = readRDS('./italk_tumor_deg_cell.rds')
th17_deg = readRDS('./italk_th17_deg_cell.rds')

cell_col = c('red', 'blue')
names(cell_col) = c('Th17', 'Tumor')

# 51 cytokin
cytosig = read.table('../BioSoftware/CytoSig/CytoSig/signature.centroid.expand',sep='\t', header = T, check.names = F)
cytosig = colnames(cytosig)
cytosig = c(cytosig,'TNFSF12')
italk_res = run_iTALK_DEG(th17_deg, tumor_deg)
par(mfrow = c(1,1), xpd = NA)
plot_data = italk_res#[1:20, ]
#plot_data = plot_data[plot_data$ligand%in%cytosig | plot_data$receptor%in%cytosig,]
wx_LRPlot_extend(plot_data,datatype='DEG',
                 cell_col=cell_col,
                 link.arr.lwd=plot_data$cell_from_logFC,
                 link.arr.width=plot_data$cell_to_logFC,
                 print.cell=TRUE
)
plot_data %>%
  ggplot(aes(x=cell_from_logFC, y=cell_to_logFC))+
    geom_point(shape=21, fill='#631879FF',color='white')+
    #geom_point(fill=ifelse(tmp2_marker$pct.1>0.25, '#631879FF', 'gray'),shape=21,color='white')+
    geom_text_repel(aes(label=ligand), max.overlaps = Inf)+
    #geom_vline(xintercept = 0.25, linetype='dashed')+
    geom_abline(slope = 1, intercept = 0)+
    theme_classic()


####
sub_srt[['cell_type']] = 'Th17'
sub_srt@meta.data[sub_srt$new_idents%in%c('Regev_tumor', 'CRC_tumor'),'cell_type'] = 'Tumor'
sub_srt[['compare_group']] = 'CRC'
sub_srt@meta.data[sub_srt$new_idents%in%c('Regev_tumor', 'Regev_Th17'),'compare_group'] = 'Regev'
# sample level
sub_srt[['sample_group']] = paste0(sub_srt$orig.ident,',', sub_srt$cell_type, ',', sub_srt$compare_group)
data = AverageExpression(sub_srt, group.by = 'sample_group')$RNA %>% log1p()
data = as.data.frame(t(data))

data$cell_type=sapply(rownames(data), function(x)strsplit(x,',')[[1]][2])
data$compare_group=sapply(rownames(data), function(x)strsplit(x, ',')[[1]][3])

deg_th17 <- DEG(data %>% filter(cell_type=='Th17'),method='Wilcox',contrast=c('CRC','Regev'))
deg_tumor<-DEG(data %>% filter(cell_type=='Tumor'),method='Wilcox',contrast=c('CRC','Regev'))

deg_th17_filter = na.omit(deg_th17[deg_th17$q.value<0.01, ])
deg_th17_filter = deg_th17_filter[!is.infinite(deg_th17_filter$logFC), ]
deg_tumor_filter = na.omit(deg_tumor[deg_tumor$q.value<0.01, ])
deg_tumor_filter = deg_tumor_filter[!is.infinite(deg_tumor_filter$logFC), ]

italk_res = run_iTALK_DEG(deg_th17_filter, deg_tumor_filter)
saveRDS(italk_res, './italk_tumor_deg_sample.rds')
plot_data = italk_res#[1:20, ]
plot_data = plot_data[plot_data$ligand%in%cytosig | plot_data$receptor%in%cytosig,]
plot_data = plot_data[order(plot_data$cell_from_logFC,decreasing=T),]
top_ligand = unique(plot_data$ligand)[1:10] 
plot_data = plot_data[plot_data$ligand%in%top_ligand, ]
saveRDS(plot_data, './italk_tumor_deg_sample_filter.rds')

par(mfrow = c(1,1), xpd = NA)
wx_LRPlot_extend(plot_data,datatype='DEG',
                 cell_col=cell_col,
                 link.arr.lwd=plot_data$cell_from_logFC,
                 link.arr.width=plot_data$cell_to_logFC,
                 print.cell=TRUE
)

tmp_res = tidyr::gather(tmp_res, key, value, -levels) %>% as.data.frame()
y = plot_data %>%
  mutate(new = paste0(ligand, receptor)) %>%
  reshape2::dcast(new ~ cell_from+comm_type, value.var = "cell_from_logFC") %>%
  replace(is.na(.), 1) %>%
  reshape2::melt()
y$a = sapply(as.vector(y$variable), function(x)strsplit(x,'_')[[1]][1])
y$b = sapply(as.vector(y$variable), function(x)strsplit(x,'_')[[1]][2])

  
x = plot_data %>%
  mutate(new = paste0(ligand, receptor)) %>%
  ggplot(aes(x=cell_from, y=new))+
  geom_point()+
  facet_wrap(~comm_type)
data = layer_data(x)
####
Tcell_srt = readRDS('./Regev_CRC_project.rds')
cell_mata = read.table('cell_meta.txt', header = T)
cluster_2_name = cell_mata$Cluster
names(cluster_2_name) = cell_mata$cluster_name
Tcell_srt[['Group']] = Tcell_srt$Cluster
Tcell_srt@meta.data[, 'Group'] = cluster_2_name[Tcell_srt$cluster_name]
Tcell_srt = subset(Tcell_srt, From=='CRC'&Group=='Tcell')
#genes = unique(c(plot_data[plot_data$cell_from=='Th17', 'ligand'], plot_data[plot_data$cell_from=='Th17', 'receptor']))
genes = unique(c(plot_data[plot_data$cell_from=='Th17', 'ligand'], plot_data[plot_data$cell_to=='Th17', 'receptor']))
a = DotPlot(Tcell_srt, features = genes,
        group.by = 'cluster_name',
        scale=F)
adata = a$data[, c('features.plot', 'id', 'pct.exp')]

# mean at sample level
library(forcats)
genes_exp = FetchData(Tcell_srt, vars = genes)
genes_exp[,c('orig.ident', 'cluster_name')] = Tcell_srt@meta.data[, c('orig.ident', 'cluster_name')]
adata2 = genes_exp %>%
  group_by(orig.ident, cluster_name) %>%
  summarise_all(list(mean)) %>%
  reshape2::melt(id=c('orig.ident', 'cluster_name')) %>%
  group_by(cluster_name, variable) %>%
  summarise(sample_mean_exp = mean(value)) %>% as.data.frame()

adata2 = merge(adata, adata2, by.x = c('features.plot', 'id'), by.y=c('variable', 'cluster_name'))

new_order = adata2[adata2$features.plot=="TNFSF12",]
new_order = new_order[order(new_order$sample_mean_exp), ]
adata2$id = factor(adata2$id, new_order$id)

adata2$features.plot= factor(adata2$features.plot, levels = c(
  'TNFSF12',  'LTA','IL22', 'IL21', 'IL1A', 'IL2',
  'TNFRSF25','ITGAV',  'ACVR1B', 'BMPR2', "BMPR1A"
))

adata2 %>%
  ggplot(aes(x=features.plot, y=id, fill=sample_mean_exp, size=pct.exp)) +
  geom_point(shape=21)+
  scale_fill_gradientn(colours = c('white','yellow', 'red'))+
  theme_classic()+
  labs(y='Cluster', x='Gene')



###
# run cytosig tumor
#receptor_marker = c('IL17RA', 'TNFRSF12A', 'LTBR', 'TNFRSF12', 'IL10RB', 'IL22RA1', 'TGFBR2', 'ITGB8')
ligand_marker = c('TWEAK',  'LTA','IL22', 'IL21', 'IL1A', 'IL2') # 
cyto_res = t(read.table('output_cell_expand_Tumor_sample_log1p.Zscore', sep='\t'))
rownames(cyto_res) = gsub('\\.', '-', rownames(cyto_res))
cyto_res = cyto_res[,ligand_marker]
ligand_marker = c('TNFSF12',  'LTA','IL22', 'IL21', 'IL1A', 'IL2')# 
colnames(cyto_res) = ligand_marker
cyto_res = melt(as.matrix(cyto_res))
colnames(cyto_res) = c('sample', 'gene',  'cytosig_activity')

th12_ligand_exp = AggregateExpression(my_th17, features = ligand_marker, group.by = 'orig.ident') $RNA %>% log1p()
th12_ligand_exp = melt(as.matrix(th12_ligand_exp))
colnames(th12_ligand_exp) = c('gene', 'sample', 'gene_exp')

ligand_cytosig_data = merge(cyto_res, th12_ligand_exp, by.x=c('gene', 'sample'), by.y=c('gene', 'sample'))

my_th17_meta = my_th17@meta.data[, c('orig.ident', 'lesions')] %>% distinct()
rownames(my_th17_meta)= my_th17_meta$orig.ident
ligand_cytosig_data$lesions = my_th17_meta[ligand_cytosig_data$sample, 'lesions']

library(ggpubr)
ligand_cytosig_data %>%
  mutate(gene=factor(gene, levels = c('TNFSF12', 'IL21','IL22',   'IL2','IL1A','LTA'))) %>%
  ggplot(aes(x=gene_exp, y=cytosig_activity))+
  geom_point(shape=21, fill="#008280FF", color='white', size=4)+
  # geom_point(shape=21, color='white', size=4)+
  facet_wrap(~gene, scales='free', nrow=2)+
  stat_cor(size=6)+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed', linewidth=1)+
  theme_bw()+
  theme(panel.grid = element_blank())

ligand_marker = c('BMP6',  'GDF11','BMP4', 'TWEAK', 'VEGFA') # 
cyto_res = t(read.table('output_cell_expand_Th17_sample.Zscore', sep='\t'))
rownames(cyto_res) = gsub('\\.', '-', rownames(cyto_res))
cyto_res = cyto_res[,ligand_marker]
ligand_marker = c('BMP6',  'GDF11','BMP4', 'TNFSF12', 'VEGFA')
colnames(cyto_res) = ligand_marker
cyto_res = melt(as.matrix(cyto_res))
colnames(cyto_res) = c('sample', 'gene',  'cytosig_activity')

Tumor_ligand_exp = AggregateExpression(subset(sub_srt, new_idents=='CRC_tumor'), features = ligand_marker, group.by = 'orig.ident') $RNA %>% log1p()
Tumor_ligand_exp = melt(as.matrix(Tumor_ligand_exp))
colnames(Tumor_ligand_exp) = c('gene', 'sample', 'gene_exp')

ligand_cytosig_data = merge(cyto_res, Tumor_ligand_exp, by.x=c('gene', 'sample'), by.y=c('gene', 'sample'))

library(ggpubr)
ligand_cytosig_data %>%
  mutate(gene=factor(gene, levels = c('BMP6', 'TNFSF12',  'BMP4', 'VEGFA','GDF11'))) %>%
  ggplot(aes(x=gene_exp, y=cytosig_activity))+
  geom_point(shape=21, fill="#008280FF", color='white', size=4)+
  # geom_point(shape=21, color='white', size=4)+
  facet_wrap(~gene, scales='free', nrow=2)+
  stat_cor(size=6)+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed', linewidth=1)+
  theme_bw()+
  theme(panel.grid = element_blank())

a = ligand_cytosig_data %>%
  filter(gene=='TNFSF12') %>%
  ggplot(aes(x=gene_exp, y=cytosig_activity))+
  geom_point(shape=21, fill="#008280FF", color='white', size=4)+
  # geom_point(shape=21, color='white', size=4)+
  labs(x='TNFSF12 exp in Tumor', y='Receptor signaling activity in Th17')+
  stat_cor(size=6)+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed', linewidth=1)+
  theme_bw()+
  theme(panel.grid = element_blank())
a
ggsave('Figure 4D_1.pdf',
       a, width=140, height=140, units='mm', dpi = 600, bg = 'transparent')


#1. mCRC vs nmCRC Th17 TWEAK expr
VlnPlot(subset(sub_srt, new_idents%in%c('CRC_Th17', 'Regev_Th17')), features = 'TNFSF12', pt.size = 0)
tmp_srt = subset(sub_srt, new_idents%in%c('CRC_Th17', 'Regev_Th17'))

tweak_exp_sample = AverageExpression(tmp_srt, features='TNFSF12', group.by = 'orig.ident')$RNA %>% log1p()
tmp_meta = tmp_srt@meta.data[,c('orig.ident', 'From')] %>% distinct()
rownames(tmp_meta) = tmp_meta$orig.ident
tmp_meta$gene_exp = tweak_exp_sample[1,tmp_meta$orig.ident]
tmp_meta$From[tmp_meta$From=='CRC'] = 'mCRC'
tmp_meta$From[tmp_meta$From=='Regev'] = 'nmCRC'


a=tmp_meta %>%
  mutate(From=factor(From, c( 'nmCRC','mCRC'))) %>%
  ggplot(aes(x=From, y=gene_exp, color=From))+
  geom_boxplot()+
  geom_jitter(width=0.1, alpha=0.5, size=2)+
  stat_compare_means(comparisons = list(c('nmCRC','mCRC')), label='p.format')+
  labs(y='TNFSF12 expression level in Th17 cells')+
  theme_bw()+
  scale_color_d3()+
    theme(axis.title.x = element_blank())
ggsave('Figure 4D_2.pdf',
       a, width=80, height=100, units='mm', dpi = 600, bg = 'transparent')

#CCR5
VlnPlot(subset(sub_srt, new_idents%in%c('CRC_Th17', 'Regev_Th17')), features = 'CCR5', pt.size = 0)
tmp_srt = subset(sub_srt, new_idents%in%c('CRC_Th17', 'Regev_Th17'))

tweak_exp_sample = AverageExpression(tmp_srt, features='CCR5', group.by = 'orig.ident')$RNA %>% log1p()
tmp_meta = tmp_srt@meta.data[,c('orig.ident', 'From')] %>% distinct()
rownames(tmp_meta) = tmp_meta$orig.ident
tmp_meta$gene_exp = tweak_exp_sample[1,tmp_meta$orig.ident]
tmp_meta$From[tmp_meta$From=='CRC'] = 'mCRC'
tmp_meta$From[tmp_meta$From=='Regev'] = 'nmCRC'


a=tmp_meta %>%
  mutate(From=factor(From, c( 'nmCRC','mCRC'))) %>%
  ggplot(aes(x=From, y=gene_exp, color=From))+
  geom_boxplot()+
  geom_jitter(width=0.1, alpha=0.5, size=2)+
  stat_compare_means(comparisons = list(c('nmCRC','mCRC')), label='p.format')+
  labs(y='CCR5 expression level in Th17 cells')+
  theme_bw()+
  scale_color_d3()+
  theme(axis.title.x = element_blank())
ggsave('Fig 8F_2.pdf',
       a, width=80, height=100, units='mm', dpi = 600, bg = 'transparent')


##Macro_c3_CD163L1的CCL4
all_srt = readRDS('Regev_CRC_project.rds')
all_srt[['new_name']] = NA
all_srt@meta.data[rownames(all_srt@meta.data[all_srt$From=='CRC', ]), 'new_name'] = all_srt@meta.data[all_srt$From=='CRC', 'cluster_name']
all_srt@meta.data[rownames(all_srt@meta.data[all_srt$From=='Regev', ]), 'new_name'] = all_srt@meta.data[all_srt$From=='Regev', 'predicted.cluster_name']
Idents(all_srt) = all_srt$new_name

tmp_srt = subset(all_srt, new_name%in%c('Macro_c3_CD163L1'))
tweak_exp_sample = AverageExpression(tmp_srt, features='CCL4', group.by = 'orig.ident')$RNA %>% log1p()
tmp_meta = tmp_srt@meta.data[,c('orig.ident', 'From')] %>% distinct()
rownames(tmp_meta) = tmp_meta$orig.ident
tmp_meta$gene_exp = tweak_exp_sample[1,tmp_meta$orig.ident]
tmp_meta$From[tmp_meta$From=='CRC'] = 'mCRC'
tmp_meta$From[tmp_meta$From=='Regev'] = 'nmCRC'


a=tmp_meta %>%
  mutate(From=factor(From, c( 'nmCRC','mCRC'))) %>%
  ggplot(aes(x=From, y=gene_exp, color=From))+
  geom_boxplot()+
  geom_jitter(width=0.1, alpha=0.5, size=2)+
  stat_compare_means(comparisons = list(c('nmCRC','mCRC')), label='p.format')+
  labs(y='CCL4 expression level in Macro_c3_CD163L1 cells')+
  theme_bw()+
  scale_color_d3()+
  theme(axis.title.x = element_blank())
ggsave('Fig 8F_1.pdf',
       a, width=80, height=100, units='mm', dpi = 600, bg = 'transparent')



###
tweak_exp_sample = AverageExpression(tmp_srt, features='TNFRSF25', group.by = 'orig.ident')$RNA %>% log1p()
tmp_meta = tmp_srt@meta.data[,c('orig.ident', 'From')] %>% distinct()
rownames(tmp_meta) = tmp_meta$orig.ident
tmp_meta$gene_exp = tweak_exp_sample[1,tmp_meta$orig.ident]
tmp_meta$From[tmp_meta$From=='CRC'] = 'mCRC'
tmp_meta$From[tmp_meta$From=='Regev'] = 'nmCRC'


VlnPlot(subset(sub_srt, new_idents%in%c('CRC_tumor', 'Regev_tumor')), features = 'TNFRSF25', pt.size = 0)
tmp_srt = subset(sub_srt, new_idents%in%c('CRC_tumor', 'Regev_tumor'))

tweak_exp_sample = AverageExpression(tmp_srt, features='TNFRSF25', group.by = 'orig.ident')$RNA %>% log1p()
tmp_meta = tmp_srt@meta.data[,c('orig.ident', 'new_idents')] %>% distinct()
rownames(tmp_meta) = tmp_meta$orig.ident
tmp_meta$gene_exp = tweak_exp_sample[1,tmp_meta$orig.ident]
tmp_meta$new_idents[tmp_meta$new_idents=='CRC_tumor'] = 'mCRC'
tmp_meta$new_idents[tmp_meta$new_idents=='Regev_tumor'] = 'nmCRC'
a=tmp_meta %>%
  mutate(new_idents=factor(new_idents, c( 'nmCRC','mCRC'))) %>%
  ggplot(aes(x=new_idents, y=gene_exp, color=new_idents))+
  geom_boxplot()+
  geom_jitter(width=0.1, alpha=0.5, size=2)+
  stat_compare_means(comparisons = list(c('nmCRC','mCRC')), label='p.format')+
  labs(y='TNFRSF25 expression level in Th17 cells')+
  theme_bw()+
  scale_color_d3()+
  theme(axis.title.x = element_blank())

EMT_score = readRDS('./EMT_score_FN14.rds')
avg_exp_fun = function(obj, features, group){
  tmp_exp = FetchData(obj, vars=features, slot='data')
  tmp_exp[, group] = obj@meta.data[, group]
  mean_exp = tmp_exp %>%
    group_by(!!sym(group)) %>%
    summarise_all(list(mean))%>% as.data.frame()
  rownames(mean_exp) = mean_exp[, 1]
  mean_exp = mean_exp[,-1, drop=F]
  return(mean_exp)
}
avg_exp_fun(My_th17_tumor, 'TNFRSF12A', 'orig.ident') 
# TWEAK expression in Th17
TNFSF12_in_TH17 = AverageExpression(my_th17, features='TNFSF12', group.by = 'orig.ident')$RNA# %>% log1p()
EMT_score$TNFSF12_exp_in_Th17 = TNFSF12_in_TH17[1,EMT_score$orig.ident]
# TNFRSF25 expression in Th17
tumor = readRDS('./tumors_filter_by_cnv_man.rds')

TNFRSF25_in_Tumor = AverageExpression(tumor, features='TNFRSF25', group.by = 'orig.ident')$RNA# %>% log1p()
EMT_score$TNFRSF25_exp_in_Tumor = TNFRSF25_in_Tumor[1,EMT_score$orig.ident]

# Th17 proportion
Tcell_srt = readRDS('./Regev_CRC_project.rds')
cell_mata = read.table('cell_meta.txt', header = T)
cluster_2_name = cell_mata$Cluster
names(cluster_2_name) = cell_mata$cluster_name
Tcell_srt[['Group']] = Tcell_srt$Cluster
Tcell_srt@meta.data[, 'Group'] = cluster_2_name[Tcell_srt$cluster_name]
Tcell_srt = subset(Tcell_srt, From=='CRC'&Group=='Tcell')

Th17_prop = Tcell_srt@meta.data %>%
  group_by(orig.ident, cluster_name) %>%
  summarise(cell_num=n()) %>%
  left_join(Tcell_srt@meta.data %>% group_by(orig.ident) %>% summarise(all_num=n())) %>%
  mutate(cell_prop=cell_num/all_num) %>% 
  as.data.frame() %>%
  filter(cluster_name=='CD4_IL17A_Th17')
rownames(Th17_prop) = Th17_prop$orig.ident
EMT_score$Th17_prop = Th17_prop[EMT_score$orig.ident, 'cell_prop']

# CD4
CD4_Tcell_srt = subset(Tcell_srt, cluster_name%in%c("CD4_CCR7_Tn-like1", "CD4_FOXP3_Treg", "CD4_IL17A_Th17", "CD4_IL7R_Tn-like2"))
CD4_Th17_prop = CD4_Tcell_srt@meta.data %>%
  group_by(orig.ident, cluster_name) %>%
  summarise(cell_num=n()) %>%
  left_join(CD4_Tcell_srt@meta.data %>% group_by(orig.ident) %>% summarise(all_num=n())) %>%
  mutate(cell_prop=cell_num/all_num) %>% 
  as.data.frame() %>%
  filter(cluster_name=='CD4_IL17A_Th17')
rownames(CD4_Th17_prop) = CD4_Th17_prop$orig.ident
EMT_score$CD4_Th17_prop = CD4_Th17_prop[EMT_score$orig.ident, 'cell_prop']

m_Tcell_srt = subset(Tcell_srt, lesions=='Liver_meta')
m_Th17_prop = m_Tcell_srt@meta.data %>%
  group_by(orig.ident, cluster_name) %>%
  summarise(cell_num=n()) %>%
  left_join(m_Tcell_srt@meta.data %>% group_by(orig.ident) %>% summarise(all_num=n())) %>%
  mutate(cell_prop=cell_num/all_num) %>% 
  as.data.frame() %>%
  filter(cluster_name=='CD4_IL17A_Th17')
rownames(m_Th17_prop) = m_Th17_prop$orig.ident
EMT_score$m_Th17_prop = m_Th17_prop[EMT_score$orig.ident, 'cell_prop']



# TNFSF12_activity_in_Tumor
cyto_res = t(read.table('output_cell_expand_Tumor_sample_log1p.Zscore', sep='\t'))
rownames(cyto_res) = gsub('\\.', '-', rownames(cyto_res))
# cyto_res = cyto_res[,ligand_marker]
EMT_score$TNFSF12_activity_in_Tumor = cyto_res[EMT_score$orig.ident, 'TWEAK']



pd = reshape2::melt(EMT_score, id=c('orig.ident', 'lesions', 'EMT_score'))

ggplot(pd, aes(x=value, y=EMT_score, color=lesions))+
    geom_point()+
    facet_wrap(~variable, scales = 'free_x')+
    geom_smooth( aes(x=value, y=EMT_score),method = 'lm', formula = 'y ~ x', se=F, inherit.aes = F)+
    stat_cor( aes(x=value, y=EMT_score), inherit.aes = F)+
  theme_bw()+
  scale_color_igv()





Tumor=readRDS("tumors_filter_by_cnv_man.rds")
Tumor_sub=Tumor[rownames(Tumor) %in% 'TNFRSF12A']
Idents(Tumor_sub) = Tumor_sub$orig.ident #TNFRSF12A
exp_by_sample = AverageExpression(Tumor_sub, assays='RNA', slot='data')[[1]]
exp_by_sample = as.data.frame(exp_by_sample)
exp_by_sample<-t(exp_by_sample)
head(exp_by_sample) 
exp_by_sample = as.data.frame(exp_by_sample)
exp_by_sample$orig.ident<-rownames(exp_by_sample)
rownames(exp_by_sample)<-NULL







######################################
Tcell_srt = readRDS('./Regev_CRC_project.rds')
cell_mata = read.table('cell_meta.txt', header = T)
cluster_2_name = cell_mata$Cluster
names(cluster_2_name) = cell_mata$cluster_name
Tcell_srt[['Group']] = Tcell_srt$Cluster
Tcell_srt@meta.data[, 'Group'] = cluster_2_name[Tcell_srt$cluster_name]

Tcell_srt = subset(Tcell_srt, From=='CRC'&Group=='Tcell')
genes = unique(c(plot_data[plot_data$cell_from=='Th17', 'ligand'], plot_data[plot_data$cell_from=='Th17', 'receptor']))
DotPlot(Tcell_srt, features = genes, group.by = 'cluster_name')+
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust=1))

# colnames(x)%in%cytosig
cytosig = read.table('../BioSoftware/CytoSig/CytoSig/signature.centroid.expand',sep='\t', header = T, check.names = F)
cytosig = colnames(cytosig)
cytosig = c(cytosig,'TNFSF12')

genes_exp = FetchData(Tcell_srt, vars = cytosig)
genes_exp[,c('orig.ident', 'cluster_name')] = Tcell_srt@meta.data[, c('orig.ident', 'cluster_name')]
library(forcats)
genes_exp %>%
  group_by(orig.ident, cluster_name) %>%
  summarise_all(list(mean)) %>%
  reshape2::melt(id=c('orig.ident', 'cluster_name')) %>%
  filter(variable=='TNFSF12') %>%
  mutate(cluster_name=fct_reorder(cluster_name, -value)) %>%
  ggplot(aes(x=cluster_name, y=value)) +
  geom_boxplot()+
  geom_jitter(width=0.2)

x = genes_exp %>%
  group_by(orig.ident, cluster_name) %>%
  summarise_all(list(mean)) %>%
  reshape2::melt(id=c('orig.ident', 'cluster_name')) %>%
  #filter(variable=='TNFSF12') %>%
  group_by(cluster_name, variable) %>%
  summarise(value=mean(value)) %>%
  filter(value!=0) %>%
  group_by(variable) %>%
  top_n(n=1, wt=value) %>%
  as.data.frame() %>%
  arrange(cluster_name)  %>%
  reshape2::dcast(cluster_name~variable)

genes_exp2 = genes_exp[, -39]
x = genes_exp2 %>%
  group_by(cluster_name) %>%
  summarise_all(list(mean)) %>% as.data.frame()
rownames(x) = x[,1]
x = x[,-1]
x[is.na(x)] = 0
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pheatmap::pheatmap(x,
                   scale = 'row',
                   #breaks = c(-2,0,2),
                   color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                             colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                   border_color='white',
                   cluster_rows = T,
                   cluster_cols = T,
                   )

cytosig = read.table('../BioSoftware/CytoSig/CytoSig/signature.centroid.expand',sep='\t', header = T, check.names = F)
cytosig = colnames(cytosig)
cytosig = c(cytosig,'TNFSF12')

cd4_name = c( "CD4_CCR7_Tn-like1","CD4_FOXP3_Treg","CD4_IL17A_Th17","CD4_IL7R_Tn-like2")
cd8_name = c("CD8_CXCL13_Tex", "CD8_GZMK_Tem", "CD8_SLC4A10_MAIT",  "Cycling_T","TYROBP_NK&NKT" )

tmp1_srt = subset(Tcell_srt, features=genes)
tmp1_srt = subset(tmp1_srt, cluster_name%in%cd4_name)
tmp1_srt[['new_ident']] = tmp1_srt$cluster_name
tmp1_srt@meta.data[tmp1_srt$new_ident!='CD4_IL17A_Th17', 'new_ident'] = 'Other_CD4'
Idents(tmp1_srt) = tmp1_srt$new_ident
tmp1_marker = FindMarkers(tmp1_srt, ident.1 = 'CD4_IL17A_Th17', ident.2 = 'Other_CD4', logfc.threshold = 0)
tmp1_marker %>% arrange(avg_log2FC)
tmp1_marker$gene = rownames(tmp1_marker)
library(ggrepel)
tmp1_marker %>%
  ggplot(aes(x=pct.1, y=pct.2, size=avg_log2FC))+
  geom_point(shape=21, fill='#631879FF',color='white')+
  geom_point(fill=ifelse(tmp1_marker$pct.1>0.25, '#631879FF', 'gray'),shape=21,color='white')+
  geom_text_repel(aes(label=ifelse(pct.1>0.25, gene, '')))+
  geom_vline(xintercept = 0.25, linetype='dashed')+
  theme_classic()
#
tmp2_srt = subset(Tcell_srt, features=genes)
tmp2_srt[['new_ident']] = tmp2_srt$cluster_name
tmp2_srt@meta.data[tmp2_srt$new_ident!='CD4_IL17A_Th17', 'new_ident'] = 'Other_CD4'
Idents(tmp2_srt) = tmp2_srt$new_ident
tmp2_marker = FindMarkers(tmp2_srt, ident.1 = 'CD4_IL17A_Th17', ident.2 = 'Other_CD4', logfc.threshold = 0)
tmp2_marker %>% arrange(avg_log2FC)
tmp2_marker$gene = rownames(tmp2_marker)
tmp2_marker = tmp2_marker %>%
  filter(gene%in%cytosig)
tmp2_marker %>%
  #filter(gene%in%cytosig) %>%
  ggplot(aes(x=pct.1, y=pct.2, size=avg_log2FC))+
  geom_point(shape=21, fill='#631879FF',color='white')+
  geom_point(fill=ifelse(tmp2_marker$pct.1>0.25, '#631879FF', 'gray'),shape=21,color='white')+
  geom_text_repel(aes(label=ifelse(pct.1>0.25, gene, '')))+
  geom_vline(xintercept = 0.25, linetype='dashed')+
  geom_abline(slope = 1, intercept = 0)+
  theme_classic()


merge_marker = tmp1_marker[,c('avg_log2FC', 'pct.1', 'pct.2')]
merge_marker[, c('avg_log2FC_all', 'pct.1_all', 'pct.2_all')] = tmp2_marker[rownames(merge_marker),c('avg_log2FC', 'pct.1', 'pct.2')]
merge_marker_FC = merge_marker[, c('avg_log2FC','avg_log2FC_all')]
colnames(merge_marker_FC) = c('Th17_vs_CD4', 'Th17_vs_OtherT')
merge_marker_FC = merge_marker_FC[rownames(merge_marker_FC)%in%cytosig, ]
merge_marker_FC = merge_marker_FC %>% arrange(-Th17_vs_OtherT, -Th17_vs_CD4,)

pheatmap::pheatmap(merge_marker_FC,
                   cluster_rows = F,cluster_cols = F)



