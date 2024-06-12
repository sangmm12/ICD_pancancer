setwd("F:/ICD/")
library(dplyr)
library(patchwork)
library(ggplotify)
library(psych)
library(stringr)
library(ggplot2)
tumor_list = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

out_put_list <- c()
for (tumor_id in tumor_list)
{
  print(tumor_id)
  dat_input <- read.table(paste0("F:/ICD/MIRNA_DATA/",tumor_id),sep='\t',check.names=F,header=T,row.names=1)
  dat_input = as.data.frame(t(dat_input))   
  na.omit(dat_input)
  dat_input$SampleName = str_sub(rownames(dat_input),end=-2)
  dat_input <- dat_input[!duplicated(dat_input$SampleName),]
  rownames(dat_input) <- make.unique(sub('\n','',dat_input$SampleName))
  dat_input <- subset(dat_input,select=-c(SampleName))
  
  
  data_all_tumor <- read.csv(paste("F:/ICD/data/ICD_tpm_exp.csv",sep=''),check.names = F)
  dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
  rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))
  dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]
  dat_one_tumor <- dat_one_tumor[grep("Tumor",dat_one_tumor$Group),]
  dat_one_tumor <- subset(dat_one_tumor,select=-c(Group,SampleName,CODE))
  one_tumor_sample <- unlist(rownames(dat_one_tumor))
  all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))
  dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]
  dat_im <- dat_input[match(all_name,rownames(dat_input)),]
  
  if(length(rownames(dat_gene))==0)
  {
    next
  }
  data.corr <- corr.test(dat_im, dat_gene, method="pearson", adjust="fdr")
  
  
  
  
  
  
  data.r = data.corr$r
  data.p = data.corr$p
  getTotal_p <- function(dc) {
    if(is.na(dc)){return(0)}
    if(dc<0.05) return(1) else return(0)
  }
  
  getTotal_r <- function(dc) {
    if(is.na(dc)){return(0)}
    if(dc<(-0.3)||dc>0.3) return(1) else return(0)
  }
  
  for(num in 10:1000)
  {
    temp.sum = as.data.frame(matrix(sapply(data.p, getTotal_p), nrow=nrow(data.p)))
    rownames(temp.sum) = rownames(data.p)
    name = names(rowSums(temp.sum)[order(-rowSums(temp.sum))])[1:num]
    temp.sum = as.data.frame(matrix(sapply(data.p, getTotal_r), nrow=nrow(data.r)))
    rownames(temp.sum) = rownames(data.r)
    name = name[which(name%in%names(rowSums(temp.sum)[order(-rowSums(temp.sum))])[1:num])]
    if(length(name)>=10)
    {
      break
    }
  }
  name = name[1:10]
  data.r = as.data.frame(data.r)[which(rownames(as.data.frame(data.r))%in%name),]
  data.p = as.data.frame(data.p)[which(rownames(as.data.frame(data.p))%in%name),]
  
  
  
  
  temp_data = data.frame()
  for(i in colnames(data.r))
  {
    for(j in rownames(data.r))
    {
      eval(parse(text=paste0("temp_row=data.frame(Row='",i,"',Col='",j,"',rvalue=","data.r[,which(colnames(data.r)=='",i,"')][which(rownames(data.r)=='",j,"')]",",pvalue=","data.p[,which(colnames(data.p)=='",i,"')][which(rownames(data.p)=='",j,"')])")))
      temp_data = rbind(temp_data,temp_row)
    }
  }
  
  pdf(paste('mirna/',tumor_id,".pdf",sep=''),width = 7,height = 10)
  p = ggplot(temp_data, aes(x = Col, y = Row, color = rvalue, size = abs(-log10(pvalue)))) +
    geom_point() +  # 使用点表示数�?
    scale_color_gradient2(low = '#623448', high = '#f8b500') +  # 调整颜色渐变,
    scale_size_continuous(range = c(1,10)) +  # 调整点的大小范围
    labs(title=tumor_id,x = "", y = "", color = "R", size = "-log10(P)") +  # 添加轴标�?
    ggtitle(tumor_id)+
    theme_bw()+
    theme(plot.title = element_text(size=30,color='black',face = "bold"),axis.text.x = element_text(angle = 45,size=20,color='black', hjust = 1),axis.text.y = element_text(size=15,color='black'))
  print(p)
  dev.off()
  
  #拼图
  eval(parse(text = paste0(tumor_id,'= p')))
  out_put_list <- append(out_put_list,tumor_id)
}
##拼图
pdf("fig1.pdf",width = length(colnames(dat_input))/40,height = length(colnames(dat_one_tumor))*length(out_put_list)/10)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))

dev.off()
