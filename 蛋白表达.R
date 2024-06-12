
library(dplyr)
library(patchwork)
library(ggplotify)
library(ggplot2)

setwd("E:/ICD/protein/")

data.corr = read.csv('protein_exp.csv',header=T,row.names=1)

data.corr[is.na(data.corr)] = 0


temp_data = data.frame()

for(i in colnames(as.data.frame(data.corr)))
{
  for(j in rownames(as.data.frame(data.corr)))
  {
    eval(parse(text=paste0("temp_row=data.frame(Row='",i,"',Col='",j,"',rvalue=","as.data.frame(data.corr)[,which(colnames(as.data.frame(data.corr))=='",i,"')][which(rownames(as.data.frame(data.corr))=='",j,"')])")))
    temp_data = rbind(temp_data,temp_row)
  }
}







pdf(paste('fig.pdf',sep=''),width = 15,height = 20)
p = ggplot(temp_data, aes(x = Col, y = Row, color = rvalue, size = rvalue)) +
  geom_point() +  # 使用点表示数�?
  scale_color_gradient(low = "white", high = "#FF9013")+
  scale_size_continuous(range = c(1,10)) +  # 调整点的大小范围
  labs(title='',x = "Gene", y = "", color = "Protein_exp", size = "Protein_exp") +  # 添加轴标�?
  ggtitle('Protein_exp')+
  theme_bw()+
  theme(plot.title = element_text(size=30,color='black',face = "bold"),axis.text.x = element_text(angle = 45,size=20,color='black', hjust = 1),axis.text.y = element_text(size=15,color='black'))
print(p)
dev.off()