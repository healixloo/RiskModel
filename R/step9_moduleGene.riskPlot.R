#install.packages("pheatmap")

library(pheatmap)
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")             
rt=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)       
rt=rt[order(rt$riskScore),]                                     

# risk curve
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file="riskScore.pdf",width = 10,height = 4)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("springgreen3",lowLength),
     rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
legend("topleft", c("High risk", "low Risk"),bty="n",pch=19,col=c("red","springgreen3"),cex=1.2)
dev.off()

# survival state plot
color=as.vector(rt$fustat)
color[color==1]="red"
color[color==0]="springgreen3"
pdf(file="survStat.pdf",width = 10,height = 4)
plot(rt$futime,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","springgreen3"),cex=1.2)
abline(v=lowLength,lty=2)
dev.off()

# risk heatmap
rt1=rt[c(3:(ncol(rt)-2))]
#rt1=log2(t(rt1)+0.001)
rt1<-t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
annotation_colors = list(
  type = c(high="red", low="springgreen3"))
rownames(annotation)=rownames(rt)
pdf(file="riskHeatmap.pdf",width = 10,height = 4)
pheatmap(rt1, 
         annotation=annotation,
          annotation_colors = annotation_colors, 
         cluster_cols = FALSE,
         fontsize_row=11,
         show_colnames = F,
         fontsize_col=3,
         scale="row",
         color = colorRampPalette(c("blue", "white", "red"))(50) )
dev.off()


