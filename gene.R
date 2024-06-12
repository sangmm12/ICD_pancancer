setwd("E:/ICD/gene_ring/")
library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
#write.csv(cyto.info,file="cyto.info.csv",quote=F) 
RCircos.Set.Core.Components(cyto.info,tracks.inside=10, tracks.outside=0 )
out.file<-"E:/ICD/gene_ring/RCircosDemoHumanGenome.pdf"#改成自己的路�?
pdf(file=out.file,height=12,width=12,compress=TRUE)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()
data(RCircos.Gene.Label.Data) 
name.col <- 4 #数据�?4�?
side <- "in" #画在基因组骨架的内侧
track.num <- 1 #基因组骨架内侧的第一个track位置上画�?
data = read.csv('E:/ICD/data/out.csv')
head(data)
RCircos.Gene.Connector.Plot(data,
                              + track.num, side)#画connector（连接基因名称和基因组位置）
track.num <- 2
RCircos.Gene.Name.Plot(data,
                         + name.col,track.num, side)#加基因名�?
#label_size=13
dev.off()