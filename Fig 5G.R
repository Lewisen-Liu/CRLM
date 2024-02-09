library(GSVA)
library(ggcorrplot)
library(stringr)
library(stringi)
library(biomaRt) 
library(dplyr)
library(igraph)
library(ggsci)

##Fig 5G#####
Genereads<-read.csv("rnaseq_genereads.csv")

rawdata=Genereads
rawdata=rawdata[-c(59310:59340,59866:59875),]
rawdata$ensembl=stri_sub(rawdata$Name,1,15)  
Ensembl_ID=rawdata$ensembl 
ensembl=useMart('ensembl',dataset = 'hsapiens_gene_ensembl') 
attributes=listAttributes(ensembl) 
ids=getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), 
          filters = 'ensembl_gene_id',
          values = Ensembl_ID, 
          mart = ensembl)
ids$hgnc_symbol[ids$hgnc_symbol==""]=NA
ids=na.omit(ids) 
colnames(ids)[1]='ensembl'
df=ids %>% left_join(rawdata,by='ensembl')

head(df)
df=df[-c(26531,6329,6345,9551,23613,20482,32980,6527),]  
rownames(df)=df$hgnc_symbol
df=df[,5:37]  

hallmark.sets = cogena::gmt2list("h.all.v7.5.1.symbols.gmt")
res.matrix = gsva(as.matrix(df), hallmark.sets, method='ssgsea', kcdf="Poisson",parallel.sz=10) 

plotdf=t(res.matrix)
head(plotdf)

plotdf=plotdf[,c('HALLMARK_G2M_CHECKPOINT',
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
                 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION')]

colnames(plotdf)=colnames(plotdf) %>% str_sub(10) 
df.cor<-cor(plotdf, method = 'pearson')
head(df.cor)
ggcorrplot(df.cor,hc.order = T,outline.color = 'gray',
           ggtheme = ggplot2::theme_dark(),
           title = "Hallmark enrichment score pearson correlation",
           tl.cex=8,tl.srt = 90) 

ggsave("Fig 5G.pdf",width = 7,height = 7,units = "in")


 



