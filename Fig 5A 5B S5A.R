library(Seurat)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(ggsci)
library(pheatmap)

malignant_tumor = readRDS("tumors_filter_by_cnv_man.rds")

###Fig 5A#####
DimPlot(malignant_tumor, group.by = 'lesions')
DimPlot(malignant_tumor, group.by = 'patient')


###Fig S5A#####
IG_name<-rownames(malignant_tumor)[grep("^IG[HLK]",rownames(malignant_tumor))]
HSP_name<-rownames(malignant_tumor)[grep("^HSP",rownames(malignant_tumor))]
keep_name<-setdiff(rownames(malignant_tumor), c(IG_name,"JCHAIN")) %>% setdiff(HSP_name)
malignant_tumor<-subset(malignant_tumor, features=keep_name)

excludeGene<-read.table("excludeGene.txt",header=TRUE)
malignant_tumor<-malignant_tumor[!rownames(malignant_tumor) %in% excludeGene$Gene]

head(malignant_tumor)
DEG_markers<-FindMarkers(malignant_tumor, ident.1="Liver_meta",group.by = "lesions",ident.2 ="Colon_tumor",logfc.threshold = 0.2) 

DEG_markers$delta<-DEG_markers$pct.1-DEG_markers$pct.2
DEG_markers$gene<-rownames(DEG_markers)
DEG_markers.filter<-DEG_markers %>% filter(delta>0) 

write.csv(DEG_markers.filter,"DEG_markers.filter.csv") #132 gene
#DEG_markers.filter=read.csv("DEG_markers.filter.csv")

malignant_tumor<-ScaleData(malignant_tumor,features = DEG_markers.filter$gene)
malignant_tumor@meta.data$CB<-rownames(malignant_tumor@meta.data) 
colanno=malignant_tumor@meta.data[,c("CB","lesions",'patient')] 
colanno=colanno%>%arrange(lesions)
rownames(colanno)=colanno$CB
colanno$CB=NULL
colanno$lesions=factor(colanno$lesions,levels = unique(colanno$lesions))


hallmark.manual=list(
  HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION=c('TIMP1','LAMC2','TGFBI','LAMA3','SDC4','PMEPA1','SPP1','ECM1','IGFBP4','CD59','FLNA','LGALS1','SERPINE2','VIM','PLOD3'),
  HALLMARK_APOPTOSIS=c('TIMP1','EMP1','GPX4','BCL2L1','DNAJA1','CLU','KRT18','GSN','ANXA1','TXNIP'),
  HALLMARK_TGF_BETA_SIGNALING=c('PMEPA1','PPP1R15A','CDKN1C'),
  HALLMARK_TNFA_SIGNALING_VIA_NFKB=c('SDC4','PMEPA1','LAMB3','PPP1R15A','FOSB','NR4A1','IER5','TUBB2A','FOSL1'),
  HALLMARK_APICAL_JUNCTION=c('LAMC2','TGFBI','LAMA3','LAMB3','MDK','MYL12B','COL17A1','CLDN4','ACTB'),
  HALLMARK_UV_RESPONSE_UP=c('DNAJA1','FOSB','NR4A1','CDKN1C','TMBIM6','DNAJB1','SELENOW'),
  HALLMARK_COAGULATION=c('TIMP1','CLU','GSN','ANXA1','CD9','APOC1'),
  HALLMARK_PROTEIN_SECRETION=c('KRT18','ARF1','TSPAN8','CD63'),
  HALLMARK_IL2_STAT5_SIGNALING=c('SPP1','ECM1','EMP1','GPX4','BCL2L1','CDKN1C','CD81','CTSZ'),
  HALLMARK_XENOBIOTIC_METABOLISM=c('IGFBP4','FBP1','TMEM176B','TMBIM6','PDK4','FABP1','APOE','SPINT2'),
  HALLMARK_KRAS_SIGNALING_UP=c('SPP1','EMP1','PPP1R15A','TMEM176B','IL2RG','DCBLD2','PCSK1N','TMEM176A'),
  HALLMARK_HYPOXIA=c('TGFBI','SDC4','PPP1R15A','CDKN1C','FBP1','TES','ANXA2','PRDX5'),
  HALLMARK_ESTROGEN_RESPONSE_LATE=c('LAMC2','IGFBP4','MDK','SLC26A2','TFF3','NBL1','KLK10','CD9'),
  HALLMARK_ESTROGEN_RESPONSE_EARLY=c('IGFBP4','KRT18','SLC26A2','TFF3','NBL1','KLK10','KRT8'),
  HALLMARK_COMPLEMENT=c('TIMP1','CD59','CLU','APOC1','CTSD','CD55','CD46'),
  HALLMARK_P53_PATHWAY=c('TXNIP','PPP1R15A','IER5','CD81','CTSD','CD82'),
  HALLMARK_MYOGENESIS=c('CLU','GSN','IFRD1','CRYAB','SOD3','TNNC2'),
  HALLMARK_MTORC1_SIGNALING=c('PPP1R15A','CD9','TES','IFRD1','CACYBP'),
  HALLMARK_INFLAMMATORY_RESPONSE=c('TIMP1','DCBLD2','CD55','CD82'),
  HALLMARK_ALLOGRAFT_REJECTION=c('TIMP1','FLNA','IL2RG','HLA-A')
)

set.seed(1)
Tumor<-AddModuleScore(object = malignant_tumor,features = hallmark.manual,name = 'hallmark')

head(Tumor)

plotdf=Tumor@meta.data[,c(16:35)]
head(plotdf)
plotdf=scale(plotdf)
plotdf<-t(plotdf)
plotdf[1:6,1:6]

rownames(plotdf)=names(hallmark.manual)%>% stringr::str_sub(10) 

plotdf=plotdf[,rownames(colanno)] 

plotdf_PT=plotdf[,colnames(plotdf)[1:16161]]
plotdf_LM=plotdf[,colnames(plotdf)[16162:50672]]
plotdf_PT=rowMeans(plotdf_PT)
plotdf_LM=rowMeans(plotdf_LM)

plotdf_df=cbind(plotdf_PT,plotdf_LM)
plotdf_df=as.data.frame(plotdf_df)
plotdf_df$diff=plotdf_df$plotdf_LM-plotdf_df$plotdf_PT
plotdf_df=plotdf_df %>% arrange(desc(diff))

plotdf=plotdf[rownames(plotdf_df),]

plotdf[plotdf>=2]=2
plotdf[plotdf < (-2)]= -2 

ann_colors=list(
  lesions=c(Colon_tumor='#1F77B4FF',Liver_meta='#FF7F0EFF'),
  patient=c(Pt01=pal_igv()(16)[1],
            Pt02=pal_igv()(16)[2],
            Pt03=pal_igv()(16)[3],
            Pt04=pal_igv()(16)[4],
            Pt05=pal_igv()(16)[5],
            Pt06=pal_igv()(16)[6],
            Pt07=pal_igv()(16)[7],
            Pt08=pal_igv()(16)[8],
            Pt09=pal_igv()(16)[9],
            Pt10=pal_igv()(16)[10],
            Pt11=pal_igv()(16)[11],
            Pt13=pal_igv()(16)[12],
            Pt15=pal_igv()(16)[13],
            Pt16=pal_igv()(16)[14],
            Pt17=pal_igv()(16)[15],
            Pt18=pal_igv()(16)[16]))
pheatmap(plotdf,cluster_rows = F,cluster_cols = F,
         show_colnames = F,
         annotation_col = colanno,
         annotation_colors = ann_colors,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         gaps_col=as.numeric(cumsum(table(colanno$lesions))),
         filename="Fig S5A.pdf",width=10,height = 7)



##fig 5B######
library(dplyr)
library(ggsci)
library(ggpubr)
library(ggplot2)

hallmark.plot<-read.csv("Hallmark_plot.csv")
hallmark.plot<-hallmark.plot[1:20,]
hallmark.plot<-hallmark.plot %>% arrange(Ratio)
hallmark.plot$GeneSet=hallmark.plot$GeneSet %>% stringr::str_sub(10) 


set.seed(1)
hallmark.plot %>% ggplot(aes(x=Ratio, y=-log10(q.value), size=Ratio, fill=-log10(q.value)))+
  geom_point(alpha=1,shape=21,color='black')+
  scale_size(range = c(4,15),name = 'Ratio')+
  scale_fill_gradient(low = '#ffffaa',high = 'red')+
  ggrepel::geom_text_repel(label=hallmark.plot$GeneSet,size=3)+
  theme_bw()
ggsave("fig 5B.pdf",width = 8,height = 7,units = "in")




