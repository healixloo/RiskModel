
#install.packages("pheatmap")

setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")                 
fdrFilter=0.05                                                   
logFCfilter=1                                                     
conNum=39                                                         
treatNum=398                                                     

rt=read.table("all.txt",sep="\t",header=T,check.names=F,row.names=1)
diffExp=read.table("diffGeneExp.txt",sep="\t",header=T,check.names=F,row.names=1)
gene=read.table("gene.txt",sep="\t",header=F)
moduleDiffAll=rt[intersect(gene[,1],row.names(rt)),]
moduleDiffGene=intersect(gene[,1],row.names(diffExp))
hmExp=diffExp[moduleDiffGene,]
moduleDiffResult=moduleDiffAll[moduleDiffGene,]

# output DEG module genes
moduleDiffResult=rbind(ID=colnames(moduleDiffResult),moduleDiffResult)
write.table(moduleDiffResult,file="moduleDiff.xls",sep="\t",col.names=F,quote=F)

# output exp of DEG module genes
moduleGeneExp=rbind(ID=colnames(hmExp),hmExp)
write.table(moduleGeneExp,file="moduleGeneExp.txt",sep="\t",col.names=F,quote=F)

# volcano
pdf(file="module_vol.pdf",height=5,width=5)
xMax=max(abs(as.numeric(as.vector(moduleDiffAll$logFC))))
yMax=max(-log10(moduleDiffAll$fdr))+1
plot(as.numeric(as.vector(moduleDiffAll$logFC)), -log10(moduleDiffAll$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(moduleDiffAll, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=0.8)
diffSub=subset(moduleDiffAll, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="blue",cex=0.8)
abline(v=0,lty=2,lwd=3)
dev.off()

# heatmap
library(pheatmap)
#hmExp=log2(hmExp+0.001)
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(diffExp)
Type=as.data.frame(Type)
pdf(file="module_heatmap.pdf",height=12,width=15)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         show_rownames = T,
         fontsize = 12,
         fontsize_row=10,
         fontsize_col=10,
         scale="row")
dev.off()


