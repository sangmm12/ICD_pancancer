rm(list=ls())
library("pheatmap")
library("jsonlite")

setwd(dir = "F:/ICD/immune_reg/receptor/cluster")
temp = list.files(pattern="*.csv")




df = read.csv(paste('F:/ICD/immune_reg/receptor/cluster/p.csv',sep=''),header=T,row.names=1)
df = replace(df,is.na(df),1)


paletteLength = 1000

getSig <- function(dc) {
  sc <- ' '
  dc = abs(dc)
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''}
  return(sc)
}

sig.mat <- matrix(sapply(as.matrix(df), getSig), nrow=nrow(as.matrix(df)))
str(sig.mat)


tran = function(temp)
{
  temp_p = c()
  for(i in temp)
  {
    if(i > 0)
    {
      if(i < 0.000000001)
      {
        i =  0.000000001
      }
      if(-log10(i)!=0){temp_p = append(temp_p,-log10(i))}else{temp_p = append(temp_p,0.99)}
    }
    else
    {
      if(i > -0.000000001)
      {
        i =  -0.000000001
      }
      temp_p = append(temp_p,log10(abs(i)))
    }
  }
  return(temp_p)
}

df_temp = sapply(df,tran)
rownames(df_temp) = rownames(df)
df = df_temp



#EST
#myColor <- colorRampPalette(c( '#00CC00', "white",'#FF0000'))(paletteLength)
#myColor <- colorRampPalette(c( '#d3f6d1', "white",'#ff5126'))(paletteLength)
#immune reg
myColor <- colorRampPalette(c( '#008B45', "white",'#C71585'))(paletteLength)
#myColor <- colorRampPalette(c( '#623448', "white",'#f8b500'))(paletteLength)
#immune check
#myColor <- colorRampPalette(c( '#D78851', "white",'#507AAF'))(paletteLength)
#myColor <- colorRampPalette(c( '#29c6cd', "white",'#CD6600'))(paletteLength)
#myColor <- colorRampPalette(c( '#29c6cd', "white",'#d01257'))(paletteLength)
# 表达
#myColor <- colorRampPalette(c( 'blue3', "white",'red3'))(paletteLength)
# #yzx
#myColor <- colorRampPalette(c( "#BEBEBE", "white","#8B6914"))(paletteLength)
#myColor <- colorRampPalette(c( "#009D91", "white","#FF6600"))(paletteLength)
# #gx
#myColor <- colorRampPalette(c( "#96CDCD", "white","#8B4789"))(paletteLength)
#myColor <- colorRampPalette(c( "#262E47", "white","#C46F0F"))(paletteLength)
# #normal tumor
# myColor <- colorRampPalette(c( "red3", "white","#0E2185"))(paletteLength)
# CIBER
#myColor <- colorRampPalette(c( '#d3f6d1', "white",'#ff5126'))(paletteLength)
#myColor <- colorRampPalette(c( '#00CC00', "white",'#FF0000'))(paletteLength)
# fenqi
#myColor <- colorRampPalette(c( '#e04462', "white",'#f08080'))(paletteLength)

myBreaks <- c(seq(min(df), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(df)/paletteLength, max(df), length.out=floor(paletteLength/2)))


pdf('pheatmap.pdf',length(colnames(df))/2,length(rownames(df))/2)

xx <- pheatmap(df,
               color=myColor,
               breaks=myBreaks,
               clustering_method="average",main='-log10(P)',number_color='black',border_color = "black", cluster_rows=F,cluster_cols=F, cellwidth = 15,cellheight = 15,display_numbers=sig.mat)
print(xx)
dev.off()

