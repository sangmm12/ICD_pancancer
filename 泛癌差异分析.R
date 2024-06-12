
setwd("F:/ICD/chayi")
cancers <- c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')


p_value <- c()
updownzong <- c()

#cancer <- c("CESC")
for(cancer in cancers){
  #读取数据
  all_gene_data <- read.csv(paste("F:/ICD/data/ZONG.csv",sep=''),header = T)
  # 去除括号和数字的函数
  remove_parentheses_and_numbers <- function(x) {
    gsub("\\(.*\\)", "", x)
  }
  
  all_gene_data$X <- remove_parentheses_and_numbers(all_gene_data$X)
  
  #选取需要的癌症数据
  subset_df <- subset(all_gene_data, X == cancer)
  
  cluster <- cbind.data.frame(subset_df$SampleName,subset_df$Group)
  
  rownames(cluster) <- cluster[,1]
  colnames(cluster) <- c("V1","V2")
  
  library(data.table)
  

  dat <- subset_df[,-c(1,2)]
  rownames(dat) <- dat$SampleName
  # 删除列名�? "SampleName" 的列
  dat <- subset(dat, select = -c(SampleName))
  dat_immu <- dat
  
  
  rownames(dat_immu) <- gsub("[\r\n]", "",rownames(dat_immu))
  
  all_name <- names(which(table(c(rownames(dat_immu),cluster[,1] ))==2))
  
  
  dat_cluster <- cluster[match(all_name,cluster[,1]),]
  
  dat_imm <- dat_immu[match(all_name,rownames(dat_immu)),]
  dat_imm <- as.data.frame(dat_imm)
  cell <- colnames(dat_imm)
  #cell <- c("B_cells_naive","B_cells_memory","Plasma_cells","T_cells_CD8","T_cells_CD4_naive","T_cells_CD4_memory_resting","T_cells_CD4_memory_activated","T_cells_follicular_helper","T_cells_regulatory_.Tregs.","T_cells_gamma_delta","NK_cells_resting","NK_cells_activated","Monocytes","Macrophages_M0","Macrophages_M1","Macrophages_M2","Dendritic_cells_resting","Dendritic_cells_activated","Mast_cells_resting","Mast_cells_activated","Eosinophils","Neutrophils")
  dat_im <- dat_imm[cell]
  
  dat_cluster <- dat_cluster[c('V2')]
  
  dat <- data.frame()
  for(coln in colnames(dat_im))
  {
    for(i in all_name)
    {
      dat <- rbind(dat,c(coln,dat_cluster[match(i,rownames(dat_cluster)),],dat_im[c(coln)][match(i,rownames(dat_im)),]))
    }
  }
  
  dat[,3] = as.numeric(dat[,3])
  
  colnames(dat) <- c("Gene","Group","value")
  
  dat$Group <- factor(dat$Group,levels=c("Tumor","Normal"))
  
  
  library(ggpubr)
  p_gene <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  if(length(p_value)==0){
    p_value <- cbind(p_value,cancer=p_gene$p)
    rownames(p_value) <- p_gene$Gene
    
  }else {
    p_value <- p_value[match(p_gene$Gene,rownames(p_value)),]
    p_value <- cbind(p_value,cancer=p_gene$p)
  }
  
  
  
  mean <- c()
  #cell<- c("CD27")
  for(cell in unique(dat$Gene)){
    
    dat1 <- subset(dat,dat$Gene==cell)
    
    dat1$Group <- factor(dat1$Group,levels=c("Tumor","Normal"))
    
    dat1 <- dat1[,-1]
    library(dplyr)
    means=group_by(dat1, Group) %>% summarise_each(funs(mean))
    means <- means[match(c("Tumor","Normal"),means$Group),]
    
    
    if(length(mean)==0){
      mean <- means
    }else {
      mean <- cbind(mean,means)
    }
    
    
    
  }
  n1 <- grep("Group",colnames(mean))
  mean1 <- mean[,-n1]
  colnames(mean1) <- unique(dat$Gene)
  rownames(mean1) <- c("Tumor","Normal")
  
  meanzong <-  t(mean1) 
  
  updown <- rep(0, times=length(rownames(meanzong)))
  
  for( i in 1:length(rownames(meanzong))){
    if( meanzong[i,1] > meanzong[i,2] ){
      updown[i]= 1
    }else{
      updown[i]= -1
    }
  }
  
  if(length(updownzong)==0){
    updownzong <- updown
  }else {
    updownzong <- data.frame(updownzong,updown)
  }
  
  
}

colnames(p_value) <- cancers

write.csv(p_value,file=paste("p_value.csv",sep=''),quote=F)


colnames(updownzong) <- cancers
rownames(updownzong) <- rownames(meanzong)

write.csv(updownzong,file= paste("updown.csv",sep=''),quote=F)








p_value <-  fread(paste("p_value.csv",sep=''))

p_value <- as.data.frame(p_value)
rownames(p_value) <- p_value[,1]
p_value <- p_value[,-1]


up_down <- updownzong[match(rownames(p_value),rownames(updownzong)),]

df <- -log10(p_value)

for(i in 1:length(rownames(df))){
  for(t in 1:length(colnames(df))){
    df[i,t] <- as.numeric(df[i,t])*as.numeric(up_down[i,t])
  }
}




library(pheatmap)


#pdf("cluster_Immunoinhibitor.pdf",width =5.5,height = 7.1)

#p <- pheatmap(-log10(p_value),  cellwidth = 15, cellheight = 15)
pdf(paste("p.pdf",sep=''),20,length(rownames(df))/1)


paletteLength = 1000
#immune
#myColor <- colorRampPalette(c("white", "#FF7C00"))(paletteLength)
#exp
#myColor <- colorRampPalette(c("white", "red"))(paletteLength)
#cell
#myColor <- colorRampPalette(c("white","blue"))(paletteLength)
#drug
#myColor <- colorRampPalette(c("white", "#660BAB"))(paletteLength)
#yzx_gx
#myColor <- colorRampPalette(c("white", "#C7007D"))(paletteLength)
max(df)
#bk <- c(seq( -max(abs(max(df)),abs(min(df))), -0.1,by=0.01),seq(0,max(abs(max(df)),abs(min(df))),by=0.01))
bk <- c(seq( -25, -0.1,by=0.01),seq(0,25,by=0.01))
myColor <- c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2))

#myBreaks <- c(seq( min(df),-min(df), length.out=floor(paletteLength/2)))

#bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))

#######################################
getSig <- function(dc) {
  sc <- ' '
  if (dc < 0.0001) {sc <- '****'}
  else if (dc < 0.001){sc <- '***'}
  else if (dc < 0.01){sc <- '**'}
  else if (dc < 0.05) {sc <- '*'}
  else{sc <- ''}
  return(sc)
}

sig.mat <- matrix(sapply(as.matrix(p_value), getSig), nrow=nrow(as.matrix(p_value)))
str(sig.mat)
########################################
xx <- pheatmap(df,
               color=myColor,
               breaks=bk,
               clustering_method="average", cluster_rows=F,cluster_cols=F, cellwidth = 30,cellheight = 30,main="-log10(p)",display_numbers=sig.mat)
#print(xx)
dev.off()


