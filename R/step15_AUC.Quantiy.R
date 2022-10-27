

#install.packages("rms")

library(rms)
library(survivalROC)
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")                         

#TCGA
riskFile="Risk.txt"
cliFile="clinical.txt"
outFile="tcga.AUC.Quantify.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli,riskScore=risk[,(ncol(risk)-1)],riskStat=risk[,(ncol(risk))])
rt$futime=rt$futime/365


## AUC extraction
outTab=data.frame()
#colnames(outTab)<-c("Time","AUC","Label")


for (i in 3:(ncol(rt)-1)) {
  for (j in 1:5) {
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i],predict.time =j, method="KM")
    outTab=rbind(outTab,cbind(Time=j,AUC=sprintf("%.3f",roc$AUC),Label=colnames(rt)[i]))
  }
}

write.table(outTab,file="tcga.AUC.Quantify.txt",sep="\t",row.names=F,quote=F)
outTab$AUC<-as.numeric(as.character(outTab$AUC))

pdf(file=outFile,height=7.5,width=9)
ggplot(outTab, aes(x=Time, y=AUC, group=Label)) +
    geom_line(aes(color=Label))+ geom_hline(yintercept=0.5, linetype="dashed",color = "grey")+xlab("Time(Year)")
dev.off()
