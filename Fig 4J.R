library(dplyr)
TCGA_data=read.delim("COADREAD_normalized_data.txt")
TCGA_clin=read.delim("COADREAD.clin.merged.picked.txt")
num=as.numeric(colnames(TCGA_data) %>% substring(14,15))
num[1]=1
group=ifelse(num%in%1:9,T,F)

TCGA_data=TCGA_data[,group]


colnames(TCGA_data)=colnames(TCGA_data) %>% substring(1,12)
colnames(TCGA_clin)=colnames(TCGA_clin) %>% toupper()

data_name=colnames(TCGA_data)
data_name=data_name[-1]
clin_name=colnames(TCGA_clin)
clin_name=clin_name[-1]
ptlist=intersect(data_name, clin_name) 

length(unique(data_name))
length(unique(clin_name))
length(unique(ptlist))


sub_TCGA_data=TCGA_data[,c(colnames(TCGA_data)[1],ptlist)]
sub_TCGA_clin=TCGA_clin[,c(colnames(TCGA_clin)[1],ptlist)]


#data extract
#SPP1 #17175
#CD44 #3355
#TNFSF12 #18438
#TNFRSF12A #18420
sub_TCGA_data_1=sub_TCGA_data[c(3355,17175,18420,18438),]


sub_TCGA_clin_1=sub_TCGA_clin[c(3,4,5),]
sub_TCGA_clin_1[is.na(sub_TCGA_clin_1)]=0


sub_TCGA_data_1_new=as.data.frame(t(sub_TCGA_data_1))
sub_TCGA_clin_1_new=as.data.frame(t(sub_TCGA_clin_1))

colnames(sub_TCGA_data_1_new)=sub_TCGA_data_1_new[1,]
sub_TCGA_data_1_new=sub_TCGA_data_1_new[-1,]
sub_TCGA_data_1_new$patient=rownames(sub_TCGA_data_1_new)
rownames(sub_TCGA_data_1_new)=NULL


colnames(sub_TCGA_clin_1_new)=sub_TCGA_clin_1_new[1,]
sub_TCGA_clin_1_new=sub_TCGA_clin_1_new[-1,]
sub_TCGA_clin_1_new$patient=rownames(sub_TCGA_clin_1_new)
rownames(sub_TCGA_clin_1_new)=NULL

sub_TCGA_clin_1_new$days_to_death=as.numeric(sub_TCGA_clin_1_new$days_to_death)
sub_TCGA_clin_1_new$days_to_last_followup=as.numeric(sub_TCGA_clin_1_new$days_to_last_followup)

sub_TCGA_clin_1_new$time=sub_TCGA_clin_1_new$days_to_death + sub_TCGA_clin_1_new$days_to_last_followup


plotdf=left_join(sub_TCGA_clin_1_new,sub_TCGA_data_1_new,by='patient')
colnames(plotdf)[6:9]=c('CD44','SPP1','TNFRSF12A','TNFSF12')


library(survival)
library(survminer)
head(plotdf)
plotdf$vital_status=as.numeric(plotdf$vital_status)
plotdf$CD44=as.numeric(plotdf$CD44)
plotdf$SPP1=as.numeric(plotdf$SPP1)
plotdf$TNFRSF12A=as.numeric(plotdf$TNFRSF12A)
plotdf$TNFSF12=as.numeric(plotdf$TNFSF12)

Surv(plotdf$time,plotdf$vital_status)


head(plotdf)
sur.cut=surv_cutpoint(plotdf,time='time',event = 'vital_status',variables = c('CD44','SPP1','TNFRSF12A','TNFSF12'))
sur.cut

boxplot(plotdf$TNFSF12)

res.cut=surv_categorize(sur.cut)
head(res.cut)

####Fig 4J####
fit=survfit(Surv(time,vital_status) ~ TNFSF12, data = res.cut)
fit
ggsurvplot(fit,
           data = res.cut,pval = T,risk.table = T, 
           palette = c('red','blue'),conf.int=F,
           ggtheme=theme_classic()) 



