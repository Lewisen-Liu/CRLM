library(dplyr)
library(DESeq2)
library(stringi)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(biomaRt) 
library(fgsea)
library(msigdbr)
library(ggplot2)
library(enrichplot)

#####Fig S5B######
Genereads<-read.csv("rnaseq_genereads.csv")
colnames(Genereads) %>% sort()

rawdata=Genereads
head(rawdata)
rawdata$CRC05LM=round((rawdata$CRC05LM1+rawdata$CRC05LM2)/2) 
rawdata$CRC06LM=round((rawdata$CRC06LM1+rawdata$CRC06LM2)/2) 
rawdata$CRC08PT=round((rawdata$CRC08PT1+rawdata$CRC08PT2)/2) 

rawdata[,c('CRC13PT', 'CRC14PT','CRC16LM','CRC18PT',  
           'CRC05LM1','CRC05LM2','CRC06LM1','CRC06LM2','CRC08PT1','CRC08PT2')] =NULL

rownames(rawdata)<-rawdata$Name
rawdata=rawdata[,3:28]

gset <- rawdata[rowMeans(rawdata)>0,] 
head(gset)
gset=gset[,colnames(gset) %>% sort()] 

group_list <- rep(c('LM','PT'),13)
group_list <- factor(group_list,levels = c("PT","LM"))

condition = group_list
subject <- factor(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13))  
coldata <- data.frame(row.names = colnames(gset), condition, subject)

dds <- DESeqDataSetFromMatrix(countData = gset,
                              colData = coldata,
                              design = ~subject +condition) 
dds$condition<- relevel(dds$condition, ref = "PT")

dds <- DESeq(dds)
nrDEG_DESeq2 <- as.data.frame(results(dds))
rld <- rlog(dds)
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),] 


Ensembl_ID=stri_sub(rownames(nrDEG_DESeq2),1,15) 
ensembl=useMart('ensembl',dataset = 'hsapiens_gene_ensembl',
                host = "https://dec2021.archive.ensembl.org") 
attributes=listAttributes(ensembl) 
ids=getBM(attributes = c('ensembl_gene_id','entrezgene_id'),
          filters = 'ensembl_gene_id',
          values = Ensembl_ID, 
          mart = ensembl)
length(unique(na.omit(ids$entrezgene_id))) 
colnames(ids)[1]='ENSEMBL'

ids$entrezgene_id[ids$entrezgene_id==""]=NA
ids=na.omit(ids) 

DESeq2_res=nrDEG_DESeq2
DESeq2_res$ENSEMBL=stri_sub(rownames(DESeq2_res),1,15)
plotdf=ids %>% left_join(DESeq2_res,by='ENSEMBL')
plotdf=plotdf %>% dplyr::arrange(dplyr::desc(log2FoldChange))

ranks=plotdf$log2FoldChange
names(ranks)=plotdf$entrezgene_id
head(ranks)
ranks=sort(ranks,decreasing = T)

geneset_hallmark=msigdbr(species = 'Homo sapiens',
                         category = 'H') %>% dplyr::select(gs_name,entrez_gene) 
head(geneset_hallmark)

set.seed(1)
Res=GSEA(ranks, TERM2GENE = geneset_hallmark,seed = T)

gseaplot(Res,geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',title = '',color='red',base_size=11,
         rel_heights=c(1,5,0.5,1),subplots=1:3,pvalue_table=T,ES_geom='line')

gseaplot2(
  Res, 
  geneSetID = 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
  title = "Epithelial Mesenchymal Transition", 
  color = "green",
  base_size = 11,
  rel_heights = c(1.5, 0.5, 0.75),
  subplots = 1:3, 
  pvalue_table = F, 
  ES_geom = "line" 
) #p=0.0016

ggsave("Fig S5B_EMT.pdf",width = 5,height = 5,units = "in")



Res2=GSEA(ranks, TERM2GENE = geneset_hallmark,seed = T,pvalueCutoff=1)

gseaplot2(
  Res2, 
  geneSetID = 'HALLMARK_APOPTOSIS',
  title = "Apoptosis", 
  color = "green",
  base_size = 11,
  rel_heights = c(1.5, 0.5, 0.75),
  subplots = 1:3, 
  pvalue_table = F, 
  ES_geom = "line" 
) 
ggsave("Fig S5B_apoptosis.pdf",width = 5,height = 5,units = "in")



gseaplot2(
  Res2, 
  geneSetID = 'HALLMARK_APICAL_JUNCTION',
  title = "Apical junction", 
  color = "green",
  base_size = 11,
  rel_heights = c(1.5, 0.5, 0.75),
  subplots = 1:3, 
  pvalue_table = F, 
  ES_geom = "line" 
)
ggsave("Fig S5B_apical.pdf",width = 5,height = 5,units = "in")


gseaplot2(
  Res2, 
  geneSetID = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
  title = "TNFa signaling via NFkB", 
  color = "green",
  base_size = 11,
  rel_heights = c(1.5, 0.5, 0.75),
  subplots = 1:3, 
  pvalue_table = F, 
  ES_geom = "line" 
)
ggsave("Fig S5B_TNFa.pdf",width = 5,height = 5,units = "in")




#####Fig S5C####

df=Res@result
library(ggpubr)
library(ggplot2)
library(ggthemes)
res<-df
res$logP<--log10(res$pvalue)
head(res)
#res$Group="not-significant"
res$Group[which(res$NES>0)]="up-regulated"
res$Group[which(res$NES<0)]="down-regulated"
table(res$Group)
res$Label<-""
res$Symbol<-rownames(res) %>% stringr::str_sub(10) 

res$Symbol<-c("Coagulation","Xenobiotic metabolism","Bile acid metabolism","Complement",                       
              "G2M","Fatty acid","E2F","EMT",
              "Peroxisome","Angiogenesis","Oxphos")

up.genes<-head(res$Symbol[which(res$Group=="up-regulated")],table(res$Group)[[2]])
down.genes<-head(res$Symbol[which(res$Group=="down-regulated")],table(res$Group)[[1]])
deg.top.genes<-c(up.genes,down.genes)
res$Label[match(deg.top.genes,res$Symbol)]<-deg.top.genes

ggscatter(res,x="NES",y="logP",
          color = "Group",
          title = "bulk RNA-seq",
          palette = c('#5050FF','#CE3D32'), 
          size = 3,
          label = res$Label,
          repel = T,
          xlab = "normalized enrichment score (NES)",
          ylab = "-log10(P value)")+theme_base()+
  geom_hline(yintercept = 1.30,linetype = "dashed")+
  geom_vline(xintercept = c(-1,1),linetype = "dashed")


ggsave("Fig S5C.pdf",width = 6,height = 5,units = "in")



