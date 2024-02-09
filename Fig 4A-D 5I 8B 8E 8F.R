library(ggplot2)
library(dplyr)
library(Seurat)
library(ggsci)
library(Seurat)
library(iTALK)
library(monocle)
library(reshape2)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggtext)
library(forcats)
library(circlize)
library(igraph)
library(forcats)
# after iTALK.R, cellphonedb/
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


#Fig 4A####
# italk
plot_data = readRDS('./italk_tumor_deg_sample_filter.rds')
par(mfrow = c(1,1), xpd = NA)
cell_col = c('#749B58FF', '#F0E685FF')
names(cell_col) = c('Th17', 'Tumor')
gene_col = rep('black', length(unique(c(plot_data$ligand,plot_data$receptor))))
names(gene_col) = unique(c(plot_data$ligand,plot_data$receptor))
gene_col['TNFSF12'] = '#BA6338FF'
gene_col['TNFRSF25'] = '#BA6338FF'
gene_col['TNFRSF12A'] = '#BA6338FF'
adata2 = read.table('./data/sample_gene_exp.txt', header = T)
##
tmp_mean_exp = c()
for(i in 1:nrow(plot_data)){
  l = plot_data[i, 'ligand']
  r = plot_data[i, 'receptor']
  f = plot_data[i, 'cell_from']
  t = plot_data[i, 'cell_to']
  if(f == 'Th17'){
    mean_exp = adata2[adata2$cluster_name=='CD4_IL17A_Th17'&adata2$variable==l, 'sample_mean_exp']
  }else{
    mean_exp = adata2[adata2$cluster_name=='CD4_IL17A_Th17'&adata2$variable==r, 'sample_mean_exp']
  }
  tmp_mean_exp = c(tmp_mean_exp, mean_exp)
}
plot_data$Th17_gene_exp = tmp_mean_exp
width_in =  150/ 25.4
height_in = 150 / 25.4
dev.off()
pdf('./fig/italk_res.pdf', width=width_in, height=height_in)
wx_LRPlot_extend(plot_data,datatype='DEG',
                 cell_col=cell_col,
                 gene_col=gene_col,
                 link.arr.lwd=plot_data$Th17_gene_exp,
                 link.arr.width=plot_data$Th17_gene_exp,
                 print.cell=TRUE
)
dev.off()




#Fig 4B#####
# cellphonedb + italk
iTALK_res = readRDS('./data/italk_tumor_deg_sample_filter.rds')
iTALK_LR = paste0(iTALK_res$ligand,'-', iTALK_res$receptor)
CRC_db = read.table('./data/mCRC_significant_means.txt', header = T, sep='\t', check.names = F)
CRC_db$LR = paste0(CRC_db$gene_a, '-', CRC_db$gene_b)
Regev_db = read.table('./data/nmCRC_significant_means.txt', header = T, sep='\t', check.names = F)
Regev_db$LR = paste0(Regev_db$gene_a, '-', Regev_db$gene_b)
Regev_db_pvalue = read.table('./data/nmCRC_pvalues.txt', header = T, sep='\t', check.names = F)
Regev_db_pvalue$LR = paste0(Regev_db_pvalue$gene_a, '-', Regev_db_pvalue$gene_b)
CRC_db_pvalue = read.table('./data/mCRC_pvalues.txt', header = T, sep='\t', check.names = F)
CRC_db_pvalue$LR = paste0(CRC_db_pvalue$gene_a, '-', CRC_db_pvalue$gene_b)
####
new_CRC_db = c()
for(i in 1:nrow(iTALK_res)){
  L = iTALK_res[i, 'ligand']
  R = iTALK_res[i, 'receptor']
  from_cell = iTALK_res[i, 'cell_from']
  to_cell = iTALK_res[i, 'cell_to']
  if(from_cell=='Th17'){
    weight = CRC_db[CRC_db$gene_a==L&CRC_db$gene_b==R, 'CD4_IL17A_Th17|Tumor']
    pvalue = CRC_db_pvalue[CRC_db_pvalue$gene_a==L&CRC_db_pvalue$gene_b==R, 'CD4_IL17A_Th17|Tumor']
  }else{
    weight = CRC_db[CRC_db$gene_a==L&CRC_db$gene_b==R, 'Tumor|CD4_IL17A_Th17']
    pvalue = CRC_db_pvalue[CRC_db_pvalue$gene_a==L&CRC_db_pvalue$gene_b==R, 'Tumor|CD4_IL17A_Th17']
  }
  #print(c(length(weight),weight, pvalue))
  if(length(weight)==0){
    weight = NA
    pvalue=NA
  }
  
  new_CRC_db = rbind(new_CRC_db, c(weight, pvalue))
}
new_CRC_db = as.data.frame(new_CRC_db)
new_CRC_db$LR = paste0(iTALK_res$ligand,'-', iTALK_res$receptor)
new_CRC_db$FT = paste0(iTALK_res$cell_from,'-', iTALK_res$cell_to)
new_CRC_db$group = 'mCRC'


new_Regev_db = c()
for(i in 1:nrow(iTALK_res)){
  L = iTALK_res[i, 'ligand']
  R = iTALK_res[i, 'receptor']
  from_cell = iTALK_res[i, 'cell_from']
  to_cell = iTALK_res[i, 'cell_to']
  if(from_cell=='Th17'){
    weight = Regev_db[Regev_db$gene_a==L&Regev_db$gene_b==R, 'cTNI05 (CD4+ IL17+)|Tumor']
    pvalue = Regev_db_pvalue[Regev_db_pvalue$gene_a==L&Regev_db_pvalue$gene_b==R, 'cTNI05 (CD4+ IL17+)|Tumor']
  }else{
    weight = Regev_db[Regev_db$gene_a==L&Regev_db$gene_b==R, 'Tumor|cTNI05 (CD4+ IL17+)']
    pvalue = Regev_db_pvalue[Regev_db_pvalue$gene_a==L&Regev_db_pvalue$gene_b==R, 'Tumor|cTNI05 (CD4+ IL17+)']
  }
  if(length(weight)==0){
    weight = NA
    pvalue=NA
  }
  #print(c(weight, pvalue))
  new_Regev_db = rbind(new_Regev_db, c(weight, pvalue))
}
new_Regev_db = as.data.frame(new_Regev_db)
new_Regev_db$LR = paste0(iTALK_res$ligand,'-', iTALK_res$receptor)
new_Regev_db$FT = paste0(iTALK_res$cell_from,'-', iTALK_res$cell_to)
new_Regev_db$group = 'nmCRC'

plot_data = as.data.frame(rbind(new_CRC_db, new_Regev_db))
colnames(plot_data) = c('Score', 'Pvalue', 'LR', 'FT', 'group')
plot_data$Pvalue[is.na(plot_data$Pvalue)]=1
plot_data$Pvalue[plot_data$Pvalue<0.01]=1e-6
plot_data$Score[is.na(plot_data$Score)] = 0
a = plot_data %>%
  # filter(!is.na(value))%>%
  mutate(LR=fct_reorder(LR, -Score, min)) %>%
  ggplot(aes(x=LR, y=group, size=-log(Pvalue), fill=Score)) +
  geom_point(shape=21)+
  scale_fill_gradientn(colours = c('yellow', 'red'))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'),
  )
a
ggsave('./fig/italk_cellphonedb.pdf',
       a, width=220, height=60, units='mm', dpi = 600, bg = 'transparent')



#Fig 4C######
# gene expression in LR
plot_data = readRDS('./italk_tumor_deg_sample_filter.rds')
genes = unique(c(plot_data[plot_data$cell_from=='Th17', 'ligand'], plot_data[plot_data$cell_to=='Th17', 'receptor']))
adata = read.table('./data/dot_genes_exp.txt', header = T)
adata2 = read.table('./data/sample_gene_exp.txt', header = T)
adata2 = merge(adata, adata2, by.x = c('features.plot', 'id'), by.y=c('variable', 'cluster_name'))
new_order = adata2[adata2$features.plot=="TNFSF12",]
new_order = new_order[order(new_order$sample_mean_exp), ]
adata2$id = factor(adata2$id, new_order$id)

adata2$features.plot= factor(adata2$features.plot, levels = c(
  'TNFSF12',  'LTA','IL22', 'IL21', 'IL1A', 'IL2',
  'TNFRSF25','ITGAV',  'ACVR1B', 'BMPR2', "BMPR1A"
))

a = adata2 %>%
  ggplot(aes(x=features.plot, y=id, fill=sample_mean_exp, size=pct.exp)) +
  geom_point(shape=21)+
  scale_fill_gradientn(colours = c('white','yellow', 'red'))+
  theme_classic()+
  labs(y='Cluster', x='Gene')

ggsave('./fig/gene_expression_in_LR.pdf', a, width=240, height=100, units='mm', dpi = 600, bg = 'transparent')



#Fig 4D####
# TNFSF12 exp in Tcell and TNFSF12 activate in tumor

ligand_marker = c('TWEAK',  'LTA','IL22', 'IL21', 'IL1A', 'IL2') # 
cyto_res = t(read.table('./data/output_cell_expand_Tumor_sample_log1p.Zscore', sep='\t'))
rownames(cyto_res) = gsub('\\.', '-', rownames(cyto_res))
cyto_res = cyto_res[,ligand_marker]
ligand_marker = c('TNFSF12',  'LTA','IL22', 'IL21', 'IL1A', 'IL2')# 
colnames(cyto_res) = ligand_marker
cyto_res = melt(as.matrix(cyto_res))
colnames(cyto_res) = c('sample', 'gene',  'cytosig_activity')
th17_ligand_exp = read.table('./data/th17_gene_exp.txt',  header = T)
ligand_cytosig_data = merge(cyto_res, th17_ligand_exp, by.x=c('gene', 'sample'), by.y=c('gene', 'sample'))

my_th17_meta = read.table('./data/th17_meta.txt',  header = T)
ligand_cytosig_data$lesions = my_th17_meta[ligand_cytosig_data$sample, 'lesions']

a = ligand_cytosig_data %>%
  filter(gene=='TNFSF12') %>%
  ggplot(aes(x=gene_exp, y=cytosig_activity))+
  geom_point(shape=21, fill="#008280FF", color='white', size=4)+
  # geom_point(shape=21, color='white', size=4)+
  labs(x='TNFSF12 exp in Th17', y='Receptor signaling activity in Tumor')+
  stat_cor(size=6)+
  geom_smooth(method = 'lm', formula = 'y ~ x', se=F, color='black', linetype='dashed', linewidth=1)+
  theme_bw()+
  theme(panel.grid = element_blank())
a
ggsave('./fig/TNFSF12_Receptor_activity.pdf',
       a, width=140, height=140, units='mm', dpi = 600, bg = 'transparent')


# TNFSF12 exp in Tumor and TNFSF12 activate in Th17
ligand_marker = c('BMP6',  'GDF11','BMP4', 'TWEAK', 'VEGFA') # 
cyto_res = t(read.table('./data/output_cell_expand_Th17_sample.Zscore', sep='\t'))
rownames(cyto_res) = gsub('\\.', '-', rownames(cyto_res))
cyto_res = cyto_res[,ligand_marker]
ligand_marker = c('BMP6',  'GDF11','BMP4', 'TNFSF12', 'VEGFA')
colnames(cyto_res) = ligand_marker
cyto_res = melt(as.matrix(cyto_res))
colnames(cyto_res) = c('sample', 'gene',  'cytosig_activity')

Tumor_ligand_exp = read.table('./data/Tumor_ligand_exp.txt', header = T)
ligand_cytosig_data = merge(cyto_res, Tumor_ligand_exp, by.x=c('gene', 'sample'), by.y=c('gene', 'sample'))

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
ggsave('./fig/TNFSF12_Receptor_activity_2.pdf',
       a, width=140, height=140, units='mm', dpi = 600, bg = 'transparent')

#Fig 8E#####
# cellphonedb
Regev_db = read.table('./data/nmCRC_T_myeloid_significant_means.txt', header = T, sep='\t', check.names = F)
Regev_db$LR = paste0(Regev_db$gene_a, '_', Regev_db$gene_b)
Regev_db_pvalue = read.table('./data/nmCRC_T_myeloid_pvalues.txt', header = T, sep='\t', check.names = F)
Regev_db_pvalue$LR = paste0(Regev_db_pvalue$gene_a, '_', Regev_db_pvalue$gene_b)

CRC_db = read.table('./data/mCRC_T_myeloid_significant_means.txt', header = T, sep='\t', check.names = F)
CRC_db$LR = paste0(CRC_db$gene_a, '_', CRC_db$gene_b)
CRC_db_pvalue = read.table('./data/mCRC_T_myeloid_pvalues.txt', header = T, sep='\t', check.names = F)
CRC_db_pvalue$LR = paste0(CRC_db_pvalue$gene_a, '_', CRC_db_pvalue$gene_b)
select_col = c('Macro_c3_CD163L1|CD4_FOXP3_Treg', 'Macro_c1_CCL18|CD4_FOXP3_Treg', 'Macro_c7_CXCL10|CD4_FOXP3_Treg',
               'Macro_c3_CD163L1|CD4_IL17A_Th17', 'Macro_c1_CCL18|CD4_IL17A_Th17', 'Macro_c7_CXCL10|CD4_IL17A_Th17'
)
CRC_db_sub = CRC_db[!((CRC_db$gene_a=='') | (CRC_db$gene_b=='')), c(select_col, 'LR')]
Regev_db_sub = Regev_db[!((Regev_db$gene_a=='') | (Regev_db$gene_b=='')), c(select_col, 'LR')]

pair_set = c('CCL4_CCR5', 'CCL3_CCR5', 'CXCL2_DPP4', 'CCL3L1_DPP4', 'CXCL12_CXCR3',
             'CXCL12_CXCR4', 'ICOSLG_ICOS', 'ADORA3_ENTPD1', 'TNF_TNFRSF1A', 'PVR_CD96')
# mCRC
CRC_db = CRC_db[match(pair_set, CRC_db$LR), select_col]
rownames(CRC_db) = pair_set
CRC_db_pvalue = CRC_db_pvalue[match(pair_set, CRC_db_pvalue$LR), select_col]
rownames(CRC_db_pvalue) = pair_set

CRC_db_plot = CRC_db %>% as.matrix() %>% reshape2::melt()
colnames(CRC_db_plot) = c('pair', 'cell', 'score')
CRC_db_pvalue_plot = CRC_db_pvalue %>% as.matrix() %>% reshape2::melt()
colnames(CRC_db_pvalue_plot) = c('pair', 'cell', 'pvalue')

CRC_db_data = merge(CRC_db_plot, CRC_db_pvalue_plot, by=c('pair', 'cell'))
CRC_db_data$pvalue[CRC_db_data$pvalue<0.05] = 'Signif'
CRC_db_data$pvalue[CRC_db_data$pvalue!='Signif'] = 'ns'
CRC_db_data$group='mCRC'
# nmCRC
Regev_db = Regev_db[match(pair_set, Regev_db$LR), select_col]
rownames(Regev_db) = pair_set
Regev_db_pvalue = Regev_db_pvalue[match(pair_set, Regev_db_pvalue$LR), select_col]
rownames(Regev_db_pvalue) = pair_set

regev_db_plot = Regev_db %>% as.matrix() %>% reshape2::melt()
colnames(regev_db_plot) = c('pair', 'cell', 'score')
Regev_db_pvalue_plot = Regev_db_pvalue %>% as.matrix() %>% reshape2::melt()
colnames(Regev_db_pvalue_plot) = c('pair', 'cell', 'pvalue')
regev_db_data = merge(regev_db_plot, Regev_db_pvalue_plot, by=c('pair', 'cell'))
regev_db_data$pvalue[regev_db_data$pvalue<0.05] = 'Signif'
regev_db_data$pvalue[regev_db_data$pvalue!='Signif'] = 'ns'
regev_db_data$group='nmCRC'



all_data = rbind(regev_db_data, CRC_db_data)
all_data$from_cell = sapply(as.vector(all_data$cell), function(x)strsplit(x,'\\|')[[1]][1])
all_data$to_cell = sapply(as.vector(all_data$cell), function(x)strsplit(x,'\\|')[[1]][2])
all_data$x = paste0(all_data$from_cell, '(',all_data$group, ')')
a=all_data %>%
  # filter(!is.na(score))%>%
  mutate(pair=factor(pair, rev(pair_set))) %>%
  mutate(from_cell = factor(from_cell, c('Macro_c3_CD163L1', 'Macro_c1_CCL18', 'Macro_c7_CXCL10'))) %>%
  ggplot(aes(x=from_cell, y=pair,  size=pvalue,fill=score)) +
  geom_point(shape=21)+
  facet_wrap(~to_cell+group, nrow = 1,)+
  scale_fill_gradientn(colours = c('yellow', 'red'))+
  scale_size_manual(values = c('Signif'=4, 'ns'=NA))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1, vjust=1),
        legend.text = element_text(size = 8),
        legend.key.height = unit(0.15, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'),
        strip.background = element_blank(),
        panel.margin = unit(0, "lines"),
        panel.border = element_rect(colour = 'gray', fill=NA, size=1)
  )

ggsave('./fig/Myeloid_T.pdf', a, width=160, height=100, units='mm', dpi = 600, bg = 'transparent')

#Fig 5I#####
# TNFRSF12A_TNFRSF25
library(ggpubr)
EMT_score = readRDS('./data/EMT_score_FN14.rds')
TNFRSF25_in_Tumor = read.table('./data/TNFRSF25exp_inTumor.txt', header = T, check.names = F)
EMT_score$TNFRSF25_exp_in_Tumor = unlist(TNFRSF25_in_Tumor[1,EMT_score$orig.ident])
pd = reshape2::melt(EMT_score, id=c('orig.ident', 'lesions', 'EMT_score'))
a = ggplot(pd, aes(x=value, y=EMT_score, color=lesions))+
  geom_point()+
  facet_wrap(~variable, scales = 'free_x')+
  geom_smooth( aes(x=value, y=EMT_score),method = 'lm', formula = 'y ~ x', se=F, inherit.aes = F)+
  stat_cor( aes(x=value, y=EMT_score), inherit.aes = F)+
  theme_bw()+
  scale_color_igv()
a
ggsave('./fig/TNFRSF12A_TNFRSF25.pdf', a, width=240, height=100, units='mm', dpi = 600, bg = 'transparent')

#Fig 8B#####
# network
cell_mata = read.table('./data/cell_meta.txt', header = T)
cluster_2_name = cell_mata$Cluster
names(cluster_2_name) = cell_mata$cluster_name
cluster_2_name = c(cluster_2_name, 'Tumor cell'='Tumor')
CRC_db = read.table('./data/Tumor_other_significant_means.txt', header = T, sep='\t', check.names = F)
CRC_db = CRC_db[, 13:ncol(CRC_db)]
net_pair = colSums(!is.na(CRC_db))
CRC_db[is.na(CRC_db)] = 0
net_weight = colSums(CRC_db)

net_data = data.frame(pair_num = net_pair, weight_sum=net_weight)
net_data$from_cell = sapply(rownames(net_data), function(x)strsplit(x, '\\|')[[1]][1])
net_data$to_cell = sapply(rownames(net_data), function(x)strsplit(x, '\\|')[[1]][2])
net_data$cell_type=cluster_2_name[net_data$from_cell]

net_data = net_data %>% left_join(net_data %>% group_by(from_cell) %>% summarise(from_cell_total_pair=sum(pair_num)))
node_info = net_data[, c('from_cell','from_cell_total_pair', 'cell_type')] %>% distinct()
rownames(node_info) = node_info$from_cell
color_node = pal_igv()(length(unique(node_info$cell_type)))
names(color_node) = unique(node_info$cell_type)
node_info$color = color_node[node_info$cell_type]
write.table(node_info, './fig/CRC_db_for_cytoscape_node.txt', sep = '\t', quote = F, row.names = F)

net_data = data.frame(pair_num = net_pair, weight_sum=net_weight)
net_data$from_cell = sapply(rownames(net_data), function(x)strsplit(x, '\\|')[[1]][1])
net_data$to_cell = sapply(rownames(net_data), function(x)strsplit(x, '\\|')[[1]][2])
net_data$cell_type=cluster_2_name[net_data$from_cell]
net_data$group = apply(net_data, 1, function(x)paste0(sort(x[3:4]), collapse = ','))
net_data = net_data[,c(1:2,5)] %>%
  group_by(group) %>%
  summarise_all(list(sum)) %>% as.data.frame()
net_data$from_cell = sapply(net_data$group, function(x)strsplit(x, ',')[[1]][1])
net_data$to_cell = sapply(net_data$group, function(x)strsplit(x, ',')[[1]][2])
net_data$weight_mean = net_data$weight_sum / net_data$pair_num
net_data = net_data %>% left_join(net_data %>% group_by(from_cell) %>% summarise(from_cell_total_pair=sum(pair_num)))

net_data$cell_type=cluster_2_name[net_data$from_cell]
net_data$to_cell_type=cluster_2_name[net_data$to_cell]

net_data = net_data[net_data$from_cell != net_data$to_cell,]
write.table(net_data, './fig//CRC_db_for_cytoscape.txt', sep = '\t', quote = F, row.names = F)


CRC_db_net2 = net_data[net_data$from_cell%in%c('Tumor cell', 'CD4_IL17A_Th17')|net_data$to_cell%in%c('Tumor cell', 'CD4_IL17A_Th17'), ]
write.table(CRC_db_net2, './fig//CRC_db_for_cytoscape_filter.txt', sep = '\t', quote = F, row.names = F)


CRC_db_net2 = net_data[net_data$from_cell%in%c('Tumor cell')|net_data$to_cell%in%c('Tumor cell'), ]
write.table(CRC_db_net2, './fig//CRC_db_for_cytoscape_filter.txt', sep = '\t', quote = F, row.names = F)

#Fig 8F#####
# CCR5 in Th17
tmp_meta = read.table('./data/ccr5_exp_in_Th17.txt', header = T)
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
a
ggsave('./fig/CCR5_in_Th17.pdf',
       a, width=80, height=100, units='mm', dpi = 600, bg = 'transparent')


# CCL4 in Macro_c3_CD163L1
tmp_meta = read.table('./data/ccl4_exp_in_Macro_c3_CD163L1.txt', header = T)
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
a
ggsave('./fig/CCL4_in_Macro_c3_CD163L1.txt.pdf',
       a, width=80, height=100, units='mm', dpi = 600, bg = 'transparent')

