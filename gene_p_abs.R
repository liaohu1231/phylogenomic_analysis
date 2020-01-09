###function of pangenome calculation
pangenome = function(x,y,z){
  for(i in 1:y){
    for (j in 1:500) {
      data2=sample(x,i,replace = T) ###从列向随机取样i列
      suma=rowSums(data2)
      sumrow=length(which(suma>0))
      write.table(sumrow,file = z,append = T,eol = '\t', row.names = F, col.names = F, quote = F)
    }
    write.table("", file = z, append = T, row.names = F, col.names = F, quote = F)
  }
}
###calculation of core-genome
coregenome=function(x,y,z){
  for(i in 1:y){
    for (j in 1:500) {
      data2=sample(x,i,replace = T) ###从列向随机取样i列
      suma=rowSums(data2)
      sumrow=length(which(suma>=i*0.99))
      write.table(sumrow,file = z,append = T,eol = '\t', row.names = F, col.names = F, quote = F)
    }
    write.table("", file = z, append = T, row.names = F, col.names = F, quote = F)
  }
}



dat=read.csv("D:/文章与资料/重分类marinobacter2/文章1015/genome/gene_presence_absence.Rtab",sep="\t",header=T,row.names=1)
a=apply(dat, 1, function(x){
  sum(x>0)/length(x)
})
###截取preudomondales
library(pheatmap)
p=pheatmap(dat[1:6000,])
order_col = p$tree_col$order
dat_name=data.frame(dat[,order_col])
dat_pse=dat_name[,57:257]
b=apply(dat_pse,1,function(x){sum(x>0)/length(x)})
dat_name=cbind(b,dat_pse)
dat_name=dat_name[order(dat_name$b,decreasing = T),]
dat_name_2=dat_name[dat_name$b>0.99,]


a=as.matrix(a)
plot(density(a[1:6000,]))

##pangenome number

#pangenome(x = dat,y = 257,z="D:/文章与资料/重分类marinobacter2/文章1015/genome/pan_accumulation_number.txt")

#绘出pan gene积累曲线
pan = read.table("D:/文章与资料/重分类marinobacter2/文章1015/genome/pan_accumulation_number.txt", sep = '\t')
pan=t(pan[,1:257])
par(mar=c(4,4,1,1),pin=c(6:5))
x=c(1:257)
boxplot(pan, col="green", outline=F,  axes=FALSE, ann=FALSE)
axis(1, at=10*1:257) #添加横坐标
axis(2,las=1,at=20000*0:200000) #添加左侧纵坐标
coregenome(x = dat,y = 257,z="D:/文章与资料/重分类marinobacter2/文章1015/genome/core_accumulation_number_2.txt")
par(new=T) #在现在的图中添加新图，不加这个命令会生成一个新图
#core = read.table("D:/文章与资料/重分类marinobacter2/文章1015/genome/core_accumulation_number.txt", sep = '\t')
#core=t(core)
boxplot(core, col = "yellow", outline=F, axes=FALSE,ann=F)
#axis(1, at=1:257) #添加横坐标
axis(4,las=0,at=500*0:6000) #添加右侧纵坐标
title(ylab = "Number of genes (pan)",line = 4)
title(xlab = "Number of genomes",line = 2)
mtext("Number of genes (core)",srt=90,side=4, line=2)
#dat2=dat[-grep(pattern="Acinetobacter",x = colnames(dat))]
#dat2=dat2[-grep(pattern="Moraxella",x = colnames(dat2))]
#dat2=dat2[-grep(pattern="Psychrobacter",x = colnames(dat2))]
#dat2=dat2[-grep(pattern="Perlucidibaca",x = colnames(dat2))]
#dat2=dat2[-grep(pattern="Kangiella",x = colnames(dat2))]
#coregenome(x = dat2,y = 203,z="D:/文章与资料/重分类marinobacter2/文章1015/genome/core_accumulation_number_203.txt")
par(new=T) #在现在的图中添加新图，不加这个命令会生成一个新图
core_203 = read.table("D:/文章与资料/重分类marinobacter2/文章1015/genome/core_accumulation_number_203.txt", sep = '\t')
core_203=t(core_203)
boxplot(core_203, col = "red", outline=F, ylab = "Number of genes (core)", axes=FALSE, ann=F)
#axis(1, at=1:257) #添加横坐标
#axis(4,las=0,at=500*0:6000) #添加左侧纵坐标
distaance=dist(t(dat[1:100000,]),method = "manhattan")
distaance=as.matrix(distaance)

pca <- function(data){
  dat <- scale(data) #标准化
  covdat <- cov(dat)  #求协方差矩阵
  eigendat <- eigen(covdat)  #求特征值、特征向量
  eigenValue <- eigendat$values  #特征值
  eigenVector <- eigendat$vectors  #特征向量
  order_value <- order(eigenValue,decreasing = T)  #由大到小排列特征值
  values <- eigenValue[order_value] 
  valueSum <- sum(values) 
  cumVar <- cumsum(values)/valueSum * 100  #计算主成分得分。
  order_vector <- eigenVector[,order_value] 
  principal <- dat %*% order_vector  #求解主成分
  return(list(PCA=principal, cumVar=cumVar))
}



library(vegan)
pcoa <- cmdscale(distaance, k = (nrow(t(dat)) - 1), eig = TRUE)
point <- data.frame(pcoa$point)
species <- wascores(pcoa$points[,1:2], t(dat))

#坐标轴解释量（前两轴）
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)

#提取样本点坐标（前两轴）
sample_site <- data.frame({pcoa$point})[1:2]
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
sample_site$group=paste(sapply(strsplit(rownames(sample_site),split = "_"),"[",1))
library(ggplot2)
p<- ggplot(sample_site, aes(PCoA1, PCoA2,group=group)) +
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +
  #geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  #geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  #geom_polygon(data = group_border, aes(fill=Group)) + #绘制多边形区域
  geom_point(aes(color = group), alpha = 0.8) + #可在这里修改点的透明度、大小
  scale_size_continuous(range = c(1,20))+
  scale_shape_manual()+  #可在这里修改点的形状
  scale_color_manual(aesthetics = "fill") + #可在这里修改点的颜色
  #scale_fill_manual(values = c('#C673FF2E', '#73D5FF2E', '#49C35A2E', '#FF985C2E')) + #可在这里修改区块的颜色
  #guides(fill = guide_legend(order = 1), shape = guide_legend(order = 2), color = guide_legend(order = 3)) #设置图例展示顺序
  
  labs(x = paste('PCoA axis1: ', round(100 * pcoa_eig[1], 2), '%'), 
       y = paste('PCoA axis2: ', round(100 * pcoa_eig[2], 2), '%'))+
  stat_ellipse(level = 0.95, show.legend = TRUE)+
  annotate('text', colour="#8766d")
p  

