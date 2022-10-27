

#install.packages("rms")

library(rms)
library(survivalROC)
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")                         

#TCGA
riskFile="Risk.txt"
cliFile="clinical.txt"
outFile="tcga.AUC.ROC.multiple.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(cli,riskScore=risk[,(ncol(risk)-1)],riskStat=risk[,(ncol(risk))])
rt$futime=rt$futime/365


## AUC extraction
outTab=data.frame()
#colnames(outTab)<-c("Time","AUC","Label","FP","TP")

for (i in 3:(ncol(rt)-1)) {
  for (j in 1:5) {
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i],predict.time =j, method="KM")
    outTab=rbind(outTab,cbind(Time=j,AUC=sprintf("%.3f",roc$AUC),Label=colnames(rt)[i],FP=roc$FP,TP=roc$TP))
  }
}

write.table(outTab,file="tcga.AUC.ROC.multiple.txt",sep="\t",row.names=F,quote=F)
outTab$AUC<-as.numeric(as.character(outTab$AUC))
outTab$FP<-as.numeric(as.character(outTab$FP))
outTab$TP<-as.numeric(as.character(outTab$TP))

pdf(file=outFile,height=7.5,width=9)
ggplot(outTab, aes(x=Time, y=AUC, group=Label)) +
    geom_line(aes(color=Label))+ geom_hline(yintercept=0.5, linetype="dashed",color = "grey")+xlab("Time(Year)")+ scale_color_brewer(palette="Set1")
dev.off()

## plot ROC
pdf(file="ROC.multiple_1year.pdf",width=8,height=8)
outTab1<-outTab[outTab$Time==1,]
outTab1$Label2<-paste(outTab1$Label,outTab1$AUC,sep="_")
ggplot(outTab1,aes(x=FP,y=TP,group=Label2,colour=factor(Label2))) +
 geom_line()+xlab("False positive rate")+ylab("True positive rate")+ scale_color_brewer(palette="Set1")+
 theme(legend.position = c(0.75, 0.2))+ labs(colour = "AUC") 
dev.off()

pdf(file="ROC.multiple_5year.pdf",width=8,height=8)
outTab5<-outTab[outTab$Time==5,]
outTab5$Label2<-paste(outTab5$Label,outTab5$AUC,sep="_")
ggplot(outTab5,aes(x=FP,y=TP,group=Label2,colour=factor(Label2))) +
 geom_line()+xlab("False positive rate")+ylab("True positive rate")+ scale_color_brewer(palette="Set1")+
 theme(legend.position = c(0.75, 0.2))+ labs(colour = "AUC") 
dev.off()


