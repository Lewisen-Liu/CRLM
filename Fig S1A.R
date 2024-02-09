library(ComplexHeatmap)
library(circlize)
library(ggsci)
#####info
a<-read.csv("CRLM_clin.csv",header = F) #including WES
#a<-read.csv("CRLM_clin_new.csv",header = F) #including Proteome
rownames(a)<-a[,1]
a<-a[-1]
b<-t(a[1,]) 
a<-a[-1,]
a<-as.matrix(a)

col_map_func <- function(dat){
  res = list()
  num = 1
  for(n in colnames(dat)){
    tmp_name = sort(unique(dat[, n]))
    tmp_col = pal_igv()(num+length(tmp_name))[num:(num+length(tmp_name))]
    names(tmp_col) = tmp_name
    res[[n]] = tmp_col
    num = num+length(tmp_name)+1
  }
  return(res)
}
col_list = col_map_func(t(a))

col_list[[2]][[1]]<-'#d8d8d8'
col_list[[3]][[1]]<-'#d8d8d8'
col_list[[4]][[1]]<-'#d8d8d8'
col_list[[5]][[1]]<-'#d8d8d8'
col_list[[6]][[1]]<-'#d8d8d8'

col_list[[1]][[3]]<-'#56aaff'
col_list[[1]][[4]]<-'#ffaaff'
col_list[[1]][[5]]<-'#ffaa56'
col_list[[1]][[6]]<-'#00bf00'

col_list[[2]][[2]]<-'#ff7b00'
col_list[[3]][[2]]<-'#007fff'
col_list[[4]][[2]]<-'#bf0000'
col_list[[5]][[2]]<-'#5cb25c'
col_list[[6]][[2]]<-'#ffcc00'

allcol = c()
for(i in rownames(a)){
  tmp_col = col_list[[i]]
  tmp_col = tmp_col[as.character(a[i,])]
  allcol = rbind(allcol, tmp_col)
}

top_anno = HeatmapAnnotation(df=b,
                             col=list(`Pt ID`=c('1'='#5050FFFF', '2'='#CE3D32FF','3'='#749B58FF','4'='#F0E685FF',
                                                '5'='#466983FF','6'='#BA6338FF','7'='#5DB1DDFF','8'='#802268FF',
                                                '9'='#6BD76BFF','10'='#D595A7FF','11'='#924822FF','12'='#837B8DFF',
                                                '13'='#C75127FF','14'='#D58F5CFF','15'='#7A65A5FF','16'='#E4AF69FF',
                                                '17'='#3B1B53FF','18'='#CDDEB7FF','19'='#612A79FF','20'='#AE1F63FF')),
                             simple_anno_size = unit(1.5,"mm"))

#Fig S1A####
Heatmap(a,
        cluster_rows = F,cluster_columns = F,
        show_column_names = 0, show_row_names = T,
        column_split = factor(b[,1]),
        column_gap = unit(1, "mm"), 
        row_split = factor(1:nrow(a)),row_title = NULL,
        width = ncol(a)*unit(2.5, "mm"), 
        height = nrow(a)*unit(5, "mm"),
        top_annotation = top_anno,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.rect(x, y, w, h, gp = gpar(fill =allcol[i,j], col = 'white'))
        },
        heatmap_legend_param = list(
          color_bar = "discrete" 
        ),
)

#Manually exclude the treated samples and retain only the non-treated ones.


