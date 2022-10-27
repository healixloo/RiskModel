
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")             
rt=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)    

deg_cluster<-rt[,c(3:(ncol(rt)-2),ncol(rt))]


## outTab strict
outTab=data.frame()
for(i in colnames(deg_cluster[,1:ncol(deg_cluster)-1])){
geneName=i
rt=rbind(expression=deg_cluster[,i],grade=deg_cluster[,ncol(deg_cluster)])
rt=as.matrix(t(rt))
wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
conGeneMeans=mean(deg_cluster[deg_cluster$risk=="low",i])
treatGeneMeans=mean(deg_cluster[deg_cluster$risk=="high",i])
logFC=log2(treatGeneMeans)-log2(conGeneMeans)
pvalue=wilcoxTest$p.value
conMed=median(deg_cluster[deg_cluster$risk=="low",i])
treatMed=median(deg_cluster[deg_cluster$risk=="high",i])
diffMed=treatMed-conMed
if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)
outTab[outTab$fdr<0.05,]
write.csv(outTab,"./DEGs_ByRisk.csv")

## outTab_2 loose
outTab_2=data.frame()
for(i in colnames(deg_cluster[,1:ncol(deg_cluster)-1])){
    geneName=i
    rt=rbind(expression=deg_cluster[,i],grade=deg_cluster[,ncol(deg_cluster)])
    rt=as.matrix(t(rt))
    wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
    conGeneMeans=mean(deg_cluster[deg_cluster$risk=="low",i])
    treatGeneMeans=mean(deg_cluster[deg_cluster$risk=="high",i])
    logFC=log2(treatGeneMeans)-log2(conGeneMeans)
    pvalue=wilcoxTest$p.value
    conMed=median(deg_cluster[deg_cluster$risk=="low",i])
    treatMed=median(deg_cluster[deg_cluster$risk=="high",i])
    diffMed=treatMed-conMed
    outTab_2=rbind(outTab_2,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
}

pValue=outTab_2[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab_2=cbind(outTab_2,fdr=fdr)
outTab_2[outTab_2$fdr<0.05,]
write.csv(outTab_2,"./DEGs_ByRisk_2.csv")