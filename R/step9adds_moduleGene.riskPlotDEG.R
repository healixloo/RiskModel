
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")

library(pheatmap)
library(dplyr)

## pheatmap
rt=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)
# clinic2<-read.table("./clinical2.txt",header = T,stringsAsFactors = F)
clinic<-read.table("./clinical.txt",header = T,stringsAsFactors = F)
rt2<-merge(rt,clinic,by.x=0,by.y = "ID",all.x = T)
pclusters<-rt2[,c(4:26,28,31:36)]
pclusters<-pclusters[order(pclusters$risk),]
anno_pclusters<-pclusters[,c(24,25,26,27,28,29,30)]
anno_pclusters<-anno_pclusters[,c(7,6,5,4,3,2,1)]
anno_pclusters$risk<-factor(anno_pclusters$risk,levels = c("low","high"))
pdf(file="riskHeatmap2.pdf",width = 10,height = 6)
pheatmap(t(pclusters[,c(1:23)]),scale="row",show_colnames=F, show_rownames=T, cluster_cols=F, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", cex=1,fontsize = 8,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,annotation_col = anno_pclusters,color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()

## DEG pvalue outTab_2 loose
## outTab_2 loose
deg_cluster<-pclusters %>% dplyr::select(-risk,risk)
deg_cluster$risk<-factor(deg_cluster$risk,levels = c("low","high"))
outTab_2=data.frame()
for(i in colnames(deg_cluster[,1:ncol(deg_cluster)-1])){
geneName=i
rt=rbind(expression=deg_cluster[,i],grade=deg_cluster[,ncol(deg_cluster)])
rt=as.matrix(t(rt))
wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
conGeneMeans=mean(deg_cluster[deg_cluster[,ncol(deg_cluster)]==levels(deg_cluster[,ncol(deg_cluster)])[1],i])
treatGeneMeans=mean(deg_cluster[deg_cluster[,ncol(deg_cluster)]==levels(deg_cluster[,ncol(deg_cluster)])[2],i])
logFC=log2(treatGeneMeans)-log2(conGeneMeans)
pvalue=wilcoxTest$p.value
conMed=median(deg_cluster[deg_cluster[,ncol(deg_cluster)]==levels(deg_cluster[,ncol(deg_cluster)])[1],i])
treatMed=median(deg_cluster[deg_cluster[,ncol(deg_cluster)]==levels(deg_cluster[,ncol(deg_cluster)])[2],i])
diffMed=treatMed-conMed
outTab_2=rbind(outTab_2,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
}
pValue=outTab_2[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab_2=cbind(outTab_2,fdr=fdr)
# outTab_2[outTab_2$fdr<0.05,]
write.csv(outTab_2,"./DEGs_ByRisk_2s.csv")

