
setwd("~/Desktop/Pro_TCGA/TCGA/data_in_2/icgc-dataset-1654693636320")
library(dplyr)
library(survival)
library("survminer")
## Load meta files and find those new samples in icgc
icgc<-read.csv("./donor.tsv",head=T,stringsAsFactors = F,sep = '\t') #tcga<-read.csv("/Users/jlu/Desktop/Pro_TCGA/TCGA/data_in/time.txt",head=T,stringsAsFactors = F,sep = '\t')
tcga_id<-read.csv("tcga.sampleID.txt",head=T,stringsAsFactors = F,sep = '\t')
tcga_id$id<-substring(tcga_id$ID,1,12)
icgc[!(icgc$submitted_donor_id %in% tcga_id$id),]
icgc_new<-icgc[!(icgc$submitted_donor_id %in% tcga_id$id),]
## Load expression matrix for icgc new samples
icgc_new_exp<-read.csv("./exp_seq.tsv",header = T, stringsAsFactors = F,sep="\t")
icgc_new_exp_new<-icgc_new_exp[substring(icgc_new_exp$submitted_sample_id,1,12) %in% icgc_new$submitted_donor_id,]
rm(icgc_new_exp)


# remove duplicated genes which are "?" and "SLC35E2"
# check duplicated genes
xxxxx<-icgc_new_exp_new[icgc_new_exp_new$submitted_sample_id=="TCGA-AA-3821-01A-01R-1022-07","gene_id"]
table(duplicated(xxxxx))
xxxxx[which(duplicated(xxxxx)==T)] # xxxxx[which(duplicated(xxxxx))] or xxxxx[duplicated(xxxxx)]
# remove duplicated genes
icgc_new_exp_new_unique<-icgc_new_exp_new[icgc_new_exp_new$gene_id!="?"&icgc_new_exp_new$gene_id!="SLC35E2",]
# check duplicated samples
icgc_new_exp_new_unique$submitted_sample_id[duplicated(icgc_new_exp_new_unique$submitted_sample_id)]
table(icgc_new_exp_new_unique$submitted_sample_id[duplicated(icgc_new_exp_new_unique$submitted_sample_id)])
yyyyy<-as.data.frame(table(icgc_new_exp_new$submitted_sample_id))
head(yyyyy[yyyyy$Freq>20531,])
dups<-as.character(yyyyy[yyyyy$Freq>20531,"Var1"])
# remove duplicated samples
icgc_new_exp_new_unique<-icgc_new_exp_new_unique[!(icgc_new_exp_new_unique$submitted_sample_id %in% dups),]
icgc_new_exp_new_unique_wide<-tidyr::spread(icgc_new_exp_new_unique[,c("gene_id","submitted_sample_id","normalized_read_count")],submitted_sample_id,normalized_read_count)
row.names(icgc_new_exp_new_unique_wide)<-icgc_new_exp_new_unique_wide$gene_id
icgc_new_exp_new_unique_wide<-icgc_new_exp_new_unique_wide[,-1]
icgc_new_exp_new_unique_wide<-apply(icgc_new_exp_new_unique_wide,2,function(x) (x/sum(x))*1000000)
colSums(icgc_new_exp_new_unique_wide[,c(2:9)])
exp_tmp <- as.data.frame(t(icgc_new_exp_new_unique_wide))
## keep only tumor samples by TCGA barcode labels
exp_tmp<-exp_tmp[substring(row.names(exp_tmp),14,15)=="01",]

# uniform the fustat and futime information by meta file
pd<-icgc
pd$fustat<-0
pd[pd$donor_vital_status=="deceased","fustat"]<-1
pd$futime<-pd$donor_interval_of_last_followup
pd[pd$fustat==1,"futime"]<-pd[pd$fustat==1,"donor_survival_time"]
pd[,c("submitted_donor_id","fustat","futime")]
pd2<-pd[,c("submitted_donor_id","fustat","futime")]

# add survival information to exp matrix
exp_tmp$sample_id<-substring(row.names(exp_tmp),1,12)
exp_tmp<-merge(exp_tmp,pd2,by.x = "sample_id", by.y = "submitted_donor_id",all.x = T)
exp_tmp$futime<-exp_tmp$futime/365
row.names(exp_tmp)<-make.names(exp_tmp$sample_id,unique=T)
row.names(exp_tmp)<-gsub("\\.","-",row.names(exp_tmp))
### Plot1 (Method1)
############ icgc
risk_model<-read.table("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in/multiCox.xls",head=T,row.names=1,sep="\t",stringsAsFactors = F)
gse_risk<-exp_tmp[,intersect(row.names(risk_model),colnames(exp_tmp))]
gse_risk[,row.names(risk_model)[!(row.names(risk_model) %in% colnames(gse_risk))]]<-0
gse_risk_2<-gse_risk
# transform IRS expression into IRS score by gene
for(i in colnames(gse_risk_2)){gse_risk_2[,i]<-gse_risk_2[,i]*risk_model[i,"coef"]}
# calculate riskScore manually
gse_risk_2$riskScore_2<-rowSums(gse_risk_2)
gse_risk_2$risk_2=as.vector(ifelse(gse_risk_2$riskScore_2>median(gse_risk_2$riskScore_2),"high","low"))
row.names(gse_risk_2)<-gsub("\\.","-",row.names(gse_risk_2))
gse_risk_2<-merge(gse_risk_2,pd2,by.x = 0, by.y = "submitted_donor_id",all.x = T)
row.names(gse_risk_2)<-gse_risk_2[,1]
gse_risk_2<-gse_risk_2[,-1]
gse_risk_2<-na.omit(gse_risk_2)
gse_risk_2$futime<-gse_risk_2$futime/365
gse_risk_2$futime<-as.numeric(gse_risk_2$futime)
gse_risk_2$fustat<-as.numeric(gse_risk_2$fustat)
diff=survdiff(Surv(futime, fustat) ~risk_2,data = gse_risk_2)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk_2, data = gse_risk_2)
ggsurvplot(fit,
data=gse_risk_2,
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

### PLOT2 (Method2)
### loading in risk model  # multi Cox model construction
rt=read.table("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in/uniSigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
rt$futime=rt$futime/365
multiCox<-readRDS("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in/multiCox.RDS")

exp_norm_tmp<-exp_tmp[,colnames(exp_tmp) %in% colnames(rt)]
modelGenes<-colnames(rt)[3:ncol(rt)]
modelGenes[!(modelGenes %in% colnames(exp_norm_tmp))]
exp_norm_tmp[,modelGenes[!(modelGenes %in% colnames(exp_norm_tmp))]]<-0
exp_norm_tmp<-exp_norm_tmp %>% select(-fustat,fustat)
exp_norm_tmp<-exp_norm_tmp %>% select(-futime,futime)
# predict risk value in new data with risk model
exp_norm_tmp$riskScore=predict(multiCox,type="risk",newdata=exp_norm_tmp)
#risk=as.vector(ifelse(riskScore>1.1893,"high","low"))
exp_norm_tmp$risk=as.vector(ifelse(exp_norm_tmp$riskScore>median(exp_norm_tmp$riskScore),"high","low"))
# generate risk file
gse_risk_norm<-exp_norm_tmp
gse_risk_norm<-gse_risk_norm[substring(row.names(gse_risk_norm),13)=="",]
## plot survival curve
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
legend.labs=c("High risk", "Low risk"),
legend.title="Risk",
xlab="Time(years)",
break.time.by = 1,
risk.table.title="",
palette=c("red", "springgreen3"),
risk.table.height=.25)
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/validation_icgc.Rhistory")
