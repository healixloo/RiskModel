#install.packages("beeswarm")

library(beeswarm)
inputFile="moduleClinical.txt"                                         
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")         
rt=read.table(inputFile,sep="\t",header=T,check.names=F)               
clinicalNum=6                                                           
pFilter=0.05                                                           

#
outTab=data.frame(gene=colnames(rt[,(clinicalNum+2):ncol(rt)]))
for(clinical in colnames(rt[,2:(clinicalNum+1)])){
    xlabel=vector()
		tab1=table(rt[,clinical])
		labelNum=length(tab1)
		dotCol=c("blue","red")
		if(labelNum==3){
			dotCol=c(2,3,4)
		}
		if(labelNum==4){
			dotCol=c(2,3,4,5)
		}
		if(labelNum>4){
			dotCol=rainbow(labelNum)
		}
		for(i in 1:labelNum){
		  xlabel=c(xlabel,names(tab1[i]) )
		}
		clinicalPvalVector=c()
		for(i in colnames(rt[,(clinicalNum+2):ncol(rt)])){
		  rt1=rbind(expression=rt[,i],clinical=rt[,clinical])
		  rt1=as.matrix(t(rt1))
		  if(labelNum==2){
		    cliTest<-t.test(expression ~ clinical, data=rt1)
		  }else{
		    cliTest<-kruskal.test(expression ~ clinical, data = rt1)}
		  pValue=cliTest$p.value
		  stat=round(cliTest$statistic,3)
		  pval=0
		  if(pValue<0.001){
			  pval=signif(pValue,4)
			  pval=format(pval, scientific = TRUE)
			}else{
			   pval=sprintf("%.03f",pValue)
		  }
		   clinicalPvalVector=c(clinicalPvalVector,paste0(stat,"(",pval,")"))
		   if(pValue<pFilter){
				   b = boxplot(expression ~ clinical, data = rt1,outline = FALSE, plot=F)
				   yMin=min(b$stats)
				   yMax = max(b$stats/5+b$stats)
				   n = ncol(b$stats)
				   outPdf=paste0(i,".",clinical,".pdf")
				   pdf(file=outPdf,width = 7,height = 5)
				   par(mar = c(4.5,6,3,3))
				   ylab=ifelse(i=="riskScore","Risk score","Gene expression")
				   boxplot(expression ~ clinical, data = rt1,names=xlabel,
				     ylab = ylab,main=paste0(i," (p=",pval,")"),xlab=clinical,
				     cex.main=1.4, cex.lab=1.4, cex.axis=1.3,ylim=c(yMin,yMax),outline = FALSE)
				   beeswarm(expression ~ clinical, data = rt1, col =dotCol, lwd=0.1,
		         pch = 16, add = TRUE, corral="wrap")
				   dev.off()
		   }
		}
		outTab=cbind(outTab,clinicalPvalVector)
}
colnames(outTab)=c("id",colnames(rt[,2:(clinicalNum+1)]))
write.table(outTab,file="clinicalCor.xls",sep="\t",row.names=F,quote=F)

