
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")
library(Biobase)
library(GEOquery)
library(stringr)
library(hgu133plus2.db)
library(survminer)
library(survival)
gse39582 <- getGEO('GSE39582',GSEMatrix=TRUE,destdir = "/Users/jlu/Desktop/Pro_TCGA/TCGA/data_in_2/gse39582")
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
### Plot2 (Method2)
### loading in risk model
rt=read.table("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in/uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
rt$futime=rt$futime/365
multiCox<-readRDS("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in/multiCox.RDS")
exp_norm<-apply(exp, 2, function(x) x/sum(x)*1000000)
exp_norm_tmp<-as.data.frame(t(exp_norm[row.names(exp_norm) %in% colnames(rt),]))
modelGenes<-colnames(rt)[3:ncol(rt)]
modelGenes[!(modelGenes %in% colnames(exp_norm_tmp))]
exp_norm_tmp[,modelGenes[!(modelGenes %in% colnames(exp_norm_tmp))]]<-0
# predict risk value in new data with risk model
exp_norm_tmp$riskScore=predict(multiCox,type="risk",newdata=exp_norm_tmp)
exp_norm_tmp$risk=as.vector(ifelse(exp_norm_tmp$riskScore>median(exp_norm_tmp$riskScore),"high","low"))
# generate risk file
table(rownames(pd)==row.names(exp_norm_tmp))
gse_risk_norm<-cbind(pd[,c("futime","fustat")],exp_norm_tmp)
## plot survival curve
## Original plot
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
#sink("fit_validation_GSE39582.txt")
#summary(fit)
#sink()
## Additional plots ## plot survival curve
pd$Stage<-str_split_fixed(pd$characteristics_ch1.4," ",2)[,2]
pd$T<-str_split_fixed(pd$characteristics_ch1.5," ",2)[,2]
pd$N<-str_split_fixed(pd$characteristics_ch1.6," ",2)[,2]
pd$M<-str_split_fixed(pd$characteristics_ch1.7," ",2)[,2]
table(rownames(pd)==row.names(exp_norm_tmp))
gse_risk_norm<-cbind(pd[,c("futime","fustat","Stage","T","N","M")],exp_norm_tmp)
gse_risk_norm$fustat<-as.numeric(gse_risk_norm$fustat)
gse_risk_norm$futime<-as.numeric(gse_risk_norm$futime)

## gse39582_Stage.pdf
diff=survdiff(Surv(futime, fustat) ~Stage,data = gse_risk_norm)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ Stage, data = gse_risk_norm)
ggsurvplot(fit,
data=gse_risk_norm,
conf.int=TRUE,
pval=paste0("p=",pValue),
pval.size=4,
risk.table=TRUE,
xlab="Time(years)",
break.time.by = 1,
risk.table.title="",
risk.table.height=.25)
sink("fit_validation_GSE39582_Stage.txt")
summary(fit)
sink()
## gse39582_T.pdf
diff=survdiff(Surv(futime, fustat) ~T,data = gse_risk_norm)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ T, data = gse_risk_norm)
ggsurvplot(fit,
data=gse_risk_norm,
conf.int=TRUE,
pval=paste0("p=",pValue),
pval.size=4,
risk.table=TRUE,
xlab="Time(years)",
break.time.by = 1,
risk.table.title="",
risk.table.height=.25)
sink("fit_validation_GSE39582_T.txt")
summary(fit)
sink()
## gse39582_N.pdf
diff=survdiff(Surv(futime, fustat) ~N,data = gse_risk_norm)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ N, data = gse_risk_norm)
ggsurvplot(fit,
data=gse_risk_norm,
conf.int=TRUE,
pval=paste0("p=",pValue),
pval.size=4,
risk.table=TRUE,
xlab="Time(years)",
break.time.by = 1,
risk.table.title="",
risk.table.height=.25)
## gse39582_M.pdf
diff=survdiff(Surv(futime, fustat) ~M,data = gse_risk_norm)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ M, data = gse_risk_norm)
ggsurvplot(fit,
data=gse_risk_norm,
conf.int=TRUE,
pval=paste0("p=",pValue),
pval.size=4,
risk.table=TRUE,
xlab="Time(years)",
break.time.by = 1,
risk.table.title="",
risk.table.height=.25)
table(gse_risk_norm$T)
table(gse_risk_norm$M)
table(pd$M)
sink("fit_validation_GSE39582_M.txt")
summary(fit)
sink()
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/Rhistory/validation_gse39582_2.Rhistory")
