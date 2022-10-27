setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")
library(Biobase)
library(GEOquery)
library(stringr)
library(hgu133plus2.db)
library(survminer)
library(survival)
gse39582 <- getGEO('GSE39582',GSEMatrix=TRUE)
exp <- exprs(gse39582[[1]])
ids=toTable(hgu133plus2SYMBOL)
exp = exp[rownames(exp) %in% ids$probe_id,]
# if one gene matched more than one probe, then choose the probe with most mean expression as the representative
tmp = by(exp,ids$symbol,function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
exp = exp[rownames(exp) %in% probes,]
ids = ids[ids$probe_id %in% probes,]
ids=ids[match(rownames(exp),ids$probe_id),]
row.names(exp)<-ids$symbol
# filtered out normal samples
pd <- pData(gse39582[[1]])
pd<-pd[pd$source_name_ch1!="Frozen tissue of non tumoral colorectal mucosa",]
pd$fustat<-str_split_fixed(pd$characteristics_ch1.13,": ",2)[,2]
pd$futime<-as.numeric(str_split_fixed(pd$characteristics_ch1.14,": ",2)[,2])/12
## remove those with fustat of N/A
pd<-pd[pd$fustat!="N/A",]
exp<-exp[,colnames(exp) %in% row.names(pd)]
head(exp[,c(1:5)])

## Plot1 (Method1)
## apply risk model
risk_model<-read.table("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in/multiCox.xls",head=T,row.names=1,sep="\t",stringsAsFactors = F)
exp_norm<-apply(exp, 2, function(x) x/sum(x)*1000000)
exp_norm_tmp<-as.data.frame(t(exp_norm[row.names(exp_norm) %in% rownames(risk_model),]))
for(i in colnames(exp_norm_tmp)){exp_norm_tmp[,i]<-exp_norm_tmp[,i]*risk_model[i,"coef"]}
exp_norm_tmp$riskScore<-rowSums(exp_norm_tmp)
exp_norm_tmp$risk=as.vector(ifelse(exp_norm_tmp$riskScore>median(exp_norm_tmp$riskScore),"high","low"))
gse_risk_norm<-cbind(pd[,c("futime","fustat")],exp_norm_tmp)
gse_risk_norm$fustat<-as.numeric(gse_risk_norm$fustat)
diff=survdiff(Surv(futime, fustat) ~risk,data = gse_risk_norm)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = gse_risk_norm)
ggsurvplot(fit,
data=gse_risk_norm,
conf.int=TRUE,
pval=paste0("p=",pValue),
pval.size=4,
risk.table=TRUE,
legend.title="Risk",
xlab="Time(years)",
break.time.by = 1,
risk.table.title="",
palette=c("red", "springgreen3"),
risk.table.height=.25)

### Plot2 (Method2)
### loading in risk model  
rt=read.table("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in/uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
rt$futime=rt$futime/365
multiCox<-readRDS("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in/multiCox.RDS")
exp_norm<-apply(exp, 2, function(x) x/sum(x)*1000000)
exp_norm_tmp<-as.data.frame(t(exp_norm[row.names(exp_norm) %in% colnames(rt),]))
modelGenes<-colnames(rt)[3:ncol(rt)]
modelGenes[!(modelGenes %in% colnames(exp_norm_tmp))]
##### 16 of 23 genes were detected in the gse dataset
exp_norm_tmp[,modelGenes[!(modelGenes %in% colnames(exp_norm_tmp))]]<-0
# predict risk value in new data with risk model
exp_norm_tmp$riskScore=predict(multiCox,type="risk",newdata=exp_norm_tmp)
exp_norm_tmp$risk=as.vector(ifelse(exp_norm_tmp$riskScore>median(exp_norm_tmp$riskScore),"high","low"))
# generate risk file
table(rownames(pd)==row.names(exp_norm_tmp))
gse_risk_norm<-cbind(pd[,c("futime","fustat")],exp_norm_tmp)
## plot survival curve
gse_risk_norm$fustat<-as.numeric(gse_risk_norm$fustat)
gse_risk_norm$futime<-as.numeric(gse_risk_norm$futime)
diff=survdiff(Surv(futime, fustat) ~risk,data = gse_risk_norm)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk, data = gse_risk_norm)
ggsurvplot(fit,
data=gse_risk_norm,
conf.int=TRUE,
pval=paste0("p=",pValue),
pval.size=4,
risk.table=TRUE,
legend.labs=c("High risk", "Low risk"),
legend.title="Risk",
xlab="Time(years)",
break.time.by = 1,
risk.table.title="",
palette=c("red", "springgreen3"),
risk.table.height=.25)
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/validation_.Rhistory")
