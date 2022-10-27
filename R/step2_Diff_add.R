
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")
library(limma)
inputFile="symbol.txt"
fdrFilter=0.05
logFCfilter=1
#conNum=39
#treatNum=398
# Read in data, Average gene when genes were duplicated in matrix
outTab=data.frame()
#grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
# filter genes with low expression
data=data[rowMeans(data)>0.2,]
data_tumor<-data[,c(40:ncol(data))]
risk=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)
grade<-as.character(risk[substring(colnames(data_tumor),1,12),"risk"])

for(i in row.names(data_tumor)){
geneName=unlist(strsplit(i,"\\|",))[1]
geneName=gsub("\\/", "_", geneName)
rt=rbind(expression=data_tumor[i,],grade=grade)
rt=as.matrix(t(rt))
wilcoxTest<-wilcox.test(as.numeric(expression) ~ grade, data=rt)
lowGeneMeans=mean(data_tumor[i,which(grade%in%"low")])
highGeneMeans=mean(data_tumor[i,which(grade%in%"high")])
logFC=log2(highGeneMeans)-log2(lowGeneMeans)
pvalue=wilcoxTest$p.value
lowMed=median(data_tumor[i,which(grade%in%"low")])
highMed=median(data_tumor[i,which(grade%in%"high")])
diffMed=highMed-lowMed
if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
outTab=rbind(outTab,cbind(gene=i,lowMean=lowGeneMeans,highMean=highGeneMeans,logFC=logFC,pValue=pvalue))
}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab=cbind(outTab,fdr=fdr)

outTab$logFC<-as.numeric(levels(outTab$logFC))[outTab$logFC]
# output
outTab$logFC[!is.finite(outTab$logFC)]<-max(outTab$logFC[is.finite(outTab$logFC)])
write.table(outTab,file="all_DEGs_byRisk.txt",sep="\t",row.names=F,quote=F)
# output table of sigDEG
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="diff_DEGs_byRisk.xls",sep="\t",row.names=F,quote=F)

# volcano
pdf(file="vol_byRisk.pdf",height=5,width=5)
xMax=max(abs(as.numeric(as.vector(outTab$logFC))))
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=0.8)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="blue",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

# heatmap input file
heatmap=rbind(ID=colnames(data_tumor[as.vector(outDiff[,1]),]),data_tumor[as.vector(outDiff[,1]),])
write.table(heatmap,file="diffGeneExp_byRisk.txt",sep="\t",col.names=F,quote=F)
# heatmap
library(pheatmap)
hmExp=data_tumor[as.vector(outDiff[,1]),]
#hmExp=log2(hmExp+0.001)
Type=grade
names(Type)=colnames(data_tumor)
Type=as.data.frame(Type)
pdf(file="heatmap_byRisk.pdf",height=12,width=15)
pheatmap(hmExp,
annotation=Type,
color = colorRampPalette(c("blue", "white", "red"))(50),
cluster_cols =F,
show_colnames = F,
show_rownames = F,
fontsize = 12,
fontsize_row=3,
fontsize_col=10,
scale="row")
dev.off()

