library(reshape2)
dat=read.csv("F:/aai/aai_summary.tsv",sep = "\t",header = 1)
dat2=dcast(dat,formula = Genome_A~Genome_B,value.var = "Mean_AAI",drop = F,fill = 0,mean)
rownames(dat2)=colnames(dat2[,2:258])
dat2=dat2[,-1]
dat2[is.na(dat2)]=0
#dat2=dat2[order(colnames(dat2),decreasing = T),order(colnames(dat2),decreasing = T)]
dat3=t(dat2)

mat=dat2+dat3

for (i in 1:257){
  if(mat[i,i]==0.00){
    mat[i,i]=mat[i,i]+100
  }
  else{
    return(mat)
  }
}

library(pheatmap)
p=pheatmap(mat,fontsize = 3,color = colorRampPalette(c("navy", "white", "firebrick3"))(80)[20:100],
         display_numbers = T)
