

#install.packages("rms")

library(rms)
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")                         

#TCGA
riskFile="Risk.txt"
cliFile="tcgaClinical.txt"
outFile="tcga.Nomogram.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
# package data
dd <- datadist(rt)
options(datadist="dd")
# generate function
rt<-rt[,colnames(rt)!="M"]
f <- cph(Surv(futime, fustat) ~ ., x=T, y=T, singular.ok=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
# build nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)), 
    lp=F, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.6, 0.4,0.2,0.1))  
# nomogram
pdf(file=outFile,height=7.5,width=11)
plot(nom)
dev.off()

# nomogram calibrate
# for one year
cal1<-calibrate(f, cmethod="KM", method="boot",u=1,m=50,B=1000)  
pdf(file="nomogram.calibration.pdf",height=7,width=7)
plot(cal1,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) +
lines(cal1[,c('mean.predicted',"KM")], type = 'b', lwd = 2, pch = 16, col = c("#2166AC")) +
mtext("")+
box(lwd = 1) +
abline(0,1,lty = 3, lwd = 2, col = c("#224444")) 
dev.off()


##################################################################
# get all nomo scores 
library(nomogramFormula)
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)),
    lp=T, funlabel=c("1-year survival", "2-year survival", "3-year survival"),
    maxscale=100,
    fun.at=c(0.99, 0.9, 0.8, 0.6, 0.4,0.2,0.1))
results <- formula_lp(nomogram = nom)
points <- points_cal(formula = results$formula, lp = f$linear.predictors)

# nomo evaluation
library(survivalROC)
cliFile2="clinical.txt"
outFile2="tcga.AUC.Quantify.pdf"
cli2=read.table(cliFile2,sep="\t",check.names=F,header=T,row.names=1)
#sameSample=intersect(row.names(cli2),row.names(risk))
cli2=cli2[sameSample,]
rt=cbind(cli2,riskScore=risk[,(ncol(risk)-1)])
rt$futime=rt$futime/365
outTab=data.frame()
rt$nomoScore<-points

outTab=data.frame()
for (i in 3:(ncol(rt))) {
  for (j in 1:5) {
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i],predict.time =j, method="KM")
    outTab=rbind(outTab,cbind(Time=j,AUC=sprintf("%.3f",roc$AUC),Label=colnames(rt)[i]))
  }
}

write.table(outTab,file="tcga.AUC.Quantify.txt",sep="\t",row.names=F,quote=F)
outTab$AUC<-as.numeric(as.character(outTab$AUC))

pdf(file=outFile2,height=7.5,width=9)
outTab$Label<-factor(outTab$Label)
ggplot(outTab, aes(x=Time, y=AUC, group=Label)) +
    geom_line(aes(color=Label))+ geom_hline(yintercept=0.5, linetype="dashed",color = "grey")+xlab("Time(Year)")+ scale_color_brewer(palette="Set1")
dev.off()



