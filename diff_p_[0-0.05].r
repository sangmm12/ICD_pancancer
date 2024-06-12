
library("pheatmap")
library("jsonlite")

setwd(dir = "E:/ICD/protein/")
temp = list.files()
myfiles = lapply(temp, read.csv)
myfiles = lapply(myfiles, na.omit)
file_nums = length(temp)
filename = sapply(strsplit(temp,"\\."),"[[",1)
for(i in filename)
{
df = read.csv(paste0(i,'.csv'),header=T,row.names=1)


pdf(paste0(i,'.pdf'),length(colnames(df))/2,length(rownames(df))/1)

myColor <- colorRampPalette(c('white','orange'))(4)
xx <- pheatmap(df,
         color=myColor,
         breaks=c(0,1,2,3),
         clustering_method="average", cluster_rows=F,cluster_cols=F, cellwidth = 20,cellheight = 20,main="protein exp")
print(xx)
dev.off()
}
