

#install.packages("rms")

library(rms)
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in/add_nomo/one_step_add_nomo2")                         

#TCGA
riskFile="Risk.txt"
cliFile="tcgaClinical.txt"
outFile="tcga.Nomogram_add.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)])
rt2<-rt
# package data
dd <- datadist(rt)
options(datadist="dd")
# generate function
rt<-rt[,colnames(rt)!="M"]
#rt<-rt[,colnames(rt)!="riskScore"] # try it with and without riskScore by masking this command or not
f <- cph(Surv(futime, fustat) ~ ., x=T, y=T, singular.ok=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)

##################################################################
# get all nomo scores including risk score
library(nomogramFormula)
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)),
    lp=T, funlabel=c("1-year survival", "2-year survival", "3-year survival"),
    maxscale=100,
    fun.at=c(0.99, 0.9, 0.8, 0.6, 0.4,0.2,0.1))
print("checkpoint1")
results <- formula_lp(nomogram = nom)
points <- points_cal(formula = results$formula, lp = f$linear.predictors)

##################################################################
# get all nomo scores not including risk score
# package data
dd <- datadist(rt2)
options(datadist="dd")
# generate function
rt2<-rt2[,colnames(rt2)!="M"]
rt2<-rt2[,colnames(rt2)!="riskScore"] # try it with and without riskScore by masking this command or not
f <- cph(Surv(futime, fustat) ~ ., x=T, y=T, singular.ok=T, surv=T, data=rt2, time.inc=1)
surv <- Survival(f)
# build nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(2, x), function(x) surv(3, x)), 
    lp=T, funlabel=c("1-year survival", "2-year survival", "3-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.6, 0.4,0.2,0.1))  
print("checkpoint2")
results <- formula_lp(nomogram = nom)
points2 <- points_cal(formula = results$formula, lp = f$linear.predictors)




# nomo evaluation
library(survivalROC)
cliFile2="clinical.txt"
outFile2="tcga.AUC.Quantify_add2.pdf"
cli2=read.table(cliFile2,sep="\t",check.names=F,header=T,row.names=1)
#sameSample=intersect(row.names(cli2),row.names(risk))
cli2=cli2[sameSample,]
rt=cbind(cli2,riskScore=risk[,(ncol(risk)-1)])
rt$futime=rt$futime/365
outTab=data.frame()
rt$nomoScore<-points
rt$nomoScore2<-points2

outTab=data.frame()
for (i in 3:(ncol(rt))) {
  for (j in 1:5) {
    roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt[,i],predict.time =j, method="KM")
    outTab=rbind(outTab,cbind(Time=j,AUC=sprintf("%.3f",roc$AUC),Label=colnames(rt)[i]))
  }
}

write.table(outTab,file="tcga.AUC.Quantify_add2.txt",sep="\t",row.names=F,quote=F)

#outTab=read.table("tcga.AUC.Quantify_add2.txt",sep="\t",check.names=F,header=T,row.names=1) ## tcga.AUC.Quantify_add2.txt contained nomoscore include and not-include risk score which is generated from codes above manually.
outTab$AUC<-as.numeric(as.character(outTab$AUC))
pdf(file=outFile2,height=7.5,width=9)
outTab$Label<-factor(outTab$Label)
ggplot(outTab, aes(x=Time, y=AUC, group=Label)) +
    geom_line(aes(color=Label))+ geom_hline(yintercept=0.5, linetype="dashed",color = "grey")+xlab("Time(Year)")+ scale_color_brewer(palette="Set1")
dev.off()



