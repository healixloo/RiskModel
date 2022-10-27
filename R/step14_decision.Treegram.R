

#install.packages("rms")

library(rms)
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")                         

#TCGA
riskFile="Risk.txt"
#cliFile="tcgaClinical.txt"
cliFile="clinical2.txt"
outFile="tcga.DecisionTree.pdf"
risk=read.table(riskFile,header=T,sep="\t",check.names=F,row.names=1)        
cli=read.table(cliFile,sep="\t",check.names=F,header=T,row.names=1)          
sameSample=intersect(row.names(cli),row.names(risk))
risk=risk[sameSample,]
cli=cli[sameSample,]
rt=cbind(futime=risk[,1],fustat=risk[,2],cli,riskScore=risk[,(ncol(risk)-1)],riskStat=risk[,(ncol(risk))])

## Decision tree
library(rpart)
library(rpart.plot)
pdf(file=outFile,height=7.5,width=11)
pfit <- rpart(Surv(futime, fustat) ~ age + gender + stage + riskStat, data = rt)
rpart.plot(pfit)
dev.off()
