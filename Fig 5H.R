library(fmsb)

#####Fig 5H#####


data <- data.frame('E'=c(5,0,3,1),
                   'H1'=c(5,0,2,0),
                   'H2'=c(5,0,1,0),
                   'H3'=c(5,0,1,0),
                   'H4'=c(5,0,0,1),
                   'M'=c(5,0,0,5))


rownames(data)[3:4] <- c('Th17->Cycling T','Th17->Treg')


# plot with default options:
radarchart(data)

# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

radarchart(data, axistype=0, 
           seg=5,
           #custom polygon
           pcol=colors_border , pfcol=colors_in , plwd=2 , plty=2,
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", cglwd=0.8, #caxislabels=seq(0,20,5)
           #custom labels
           vlcex=1.5     #, centerzero=T
)
# Add a legend
legend(x=0.7, y=1.3, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=c(rgb(0.2,0.5,0.5,0.8), rgb(0.8,0.2,0.5,0.8)) , text.col = "black", cex=1.2, pt.cex=3)

####manually saved, width = 7,height = 5 




