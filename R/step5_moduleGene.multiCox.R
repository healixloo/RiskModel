#install.packages('survival')

library(survival)                                        
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")       
rt=read.table("uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)   
rt$futime=rt$futime/365

# multi Cox model construction
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
#multiCox=step(multiCox,direction = "both") # remove duplicated genes with high ocrrelation
multiCoxSum=summary(multiCox)
saveRDS(multiCox,"multiCox.RDS")
# output model parametre
outTab=data.frame()
outTab=cbind(
             coef=multiCoxSum$coefficients[,"coef"],
             HR=multiCoxSum$conf.int[,"exp(coef)"],
             HR.95L=multiCoxSum$conf.int[,"lower .95"],
             HR.95H=multiCoxSum$conf.int[,"upper .95"],
             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

# output risk value of patients
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
    file="risk.txt",
    sep="\t",
    quote=F,
    row.names=F)


