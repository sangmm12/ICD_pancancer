setwd("F:/ICD/")
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplotify)
library(psych)
library(ggside)
library(ggpubr)


tumor_list = c('ACC',"ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

other_data = read.csv("F:/ICD/lasso/risk.csv",row.names=1,header=T,check.names=F)

out_put_list <- c()
for (tumor_id in tumor_list)
{
  print(tumor_id)
  
  dat_input <- read.csv(paste0("drug/data/",tumor_id,".csv"),check.names = T,row.names=1,header=T)
  dat_input <- dat_input[!duplicated(rownames(dat_input)),]
  
  
  
  #显著�?
  dat <- data.frame(check.names = F)
  other_file <- other_data[grep(tumor_id,other_data$CODE),]
  other_file <- subset(other_file,select=-c(CODE))
  allname <- names(which(table(c(rownames(other_file),rownames(dat_input)))==2))
  if(length(allname)==0){next}
  for(gene_name in colnames(dat_input))
  {
    for(i in allname)
    {
      dat <- rbind(dat,c(gene_name,other_file[match(i,rownames(other_file)),],dat_input[c(gene_name)][match(i,rownames(dat_input)),]))
    }
    
  }
  colnames(dat) <- c("Gene","Group","value")
  dat[,3] = as.numeric(dat[,3])
  dat <- na.omit(dat)
  xx <- compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
  
  
  
  #线�?
  data_all_tumor <- read.csv(paste("F:/ICD/data/ICD_tpm_exp.csv",sep=''),check.names = F)
  dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
  if (length(colnames(dat_one_tumor))==0)
  {
    next
  }
  rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))
  dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]
  dat_one_tumor <- dat_one_tumor[grep("Tumor",dat_one_tumor$Group),]
  dat_one_tumor <- subset(dat_one_tumor,select=-c(Group,SampleName,CODE))
  one_tumor_sample <- unlist(rownames(dat_one_tumor))
  all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))
  dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]
  dat_im <- dat_input[match(all_name,rownames(dat_input)),]
  data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")
  temp_data = data.frame()
  for(i in colnames(as.data.frame(data.corr$r)))
  {
    for(j in rownames(as.data.frame(data.corr$r)))
    {
      eval(parse(text=paste0("temp_row=data.frame(Row='",i,"',Col='",j,"',rvalue=","as.data.frame(data.corr$r)[,which(colnames(as.data.frame(data.corr$r))=='",i,"')][which(rownames(as.data.frame(data.corr$r))=='",j,"')]",",pvalue=","as.data.frame(data.corr$p)[,which(colnames(as.data.frame(data.corr$p))=='",i,"')][which(rownames(as.data.frame(data.corr$p))=='",j,"')])")))
      temp_data = rbind(temp_data,temp_row)
    }
  }
  
  # convert_p = function(x,y)
  # {
  #   temp_c = c()
  #   for(i in x)
  #   {
  #     for(j in 1:y)
  #     {
  #       temp_c = append(temp_c,i)
  #     }
  #   }
  #   return(temp_c)
  # }
  # 
  # temp_data$Label = convert_p(xx$p.signif,length(colnames(dat_one_tumor)))
  
  
  
  
  
  #geom_text(aes(label = Label), hjust = 0, vjust = 0, size = 5)
  pdf(paste('drug/corr/',tumor_id,".pdf",sep=''),width = 15,height = 5)
  p = ggplot(temp_data, aes(x = Col, y = Row, color = rvalue, size = abs(-log10(pvalue)))) +
    geom_point() +  # 使用点表示数�?
    scale_color_gradient2(low = "blue3", high = "red3") +  # 调整颜色渐变
    scale_size_continuous(range = c(1,10)) +  # 调整点的大小范围
    labs(title=tumor_id,x = "", y = "", color = "R", size = "-log10(P)") +  # 添加轴标�?
    ggtitle(tumor_id)+
    geom_ysidelabel(aes(label=p,x=x,y=y),data=data.frame(p=xx$p.signif[order(xx$Gene)],x=rep(xx$p.signif[order(xx$Gene)],length(xx$p.signif)),y=seq(1:length(xx$p.signif))),inherit.aes=F,label.size=1)+
    theme_bw()+
    theme(plot.title = element_text(size=30,color='black',face = "bold"),axis.text.x = element_text(angle = 45,size=20,color='black', hjust = 1),axis.text.y = element_text(size=15,color='black'))
  print(p)
  dev.off()
  
  #拼图
  eval(parse(text = paste0(tumor_id,'= p')))
  out_put_list <- append(out_put_list,tumor_id)
}
##拼图
pdf("fig1.pdf",width = length(colnames(dat_input))*6,height = length(colnames(dat_one_tumor))*length(out_put_list)/25)

eval(parse(text = paste('print(cowplot::plot_grid(',paste(out_put_list, collapse = ","),', ncol=4))',sep='')))

dev.off()
