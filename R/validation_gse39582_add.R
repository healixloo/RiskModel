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
exp_norm<-apply(exp, 2, function(x) x/sum(x)*1000000)
tgenes<-read.csv("./all.comparasion.genes",header = F,stringsAsFactors = F)
exp_norm_tmp<-as.data.frame(t(exp_norm[row.names(exp_norm) %in% tgenes$V1,]))
gse_risk_norm<-cbind(pd[,c("futime","fustat")],exp_norm_tmp)
gse_risk_norm$fustat<-as.numeric(gse_risk_norm$fustat)

SurvivalByGene<-function(gse_risk_norm,tgene){
gse_risk_norm$gene=as.vector(ifelse(gse_risk_norm[,tgene]>median(gse_risk_norm[,tgene]),"high","low"))
diff=survdiff(Surv(futime, fustat) ~gene,data = gse_risk_norm)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ gene, data = gse_risk_norm)
p<-ggsurvplot(fit,
data=gse_risk_norm,
conf.int=TRUE,
pval=paste0("p=",pValue),
pval.size=4,
risk.table=TRUE,
legend.title=tgene,
xlab="Time(years)",
break.time.by = 1,
risk.table.title="",
palette=c("red", "springgreen3"),
risk.table.height=.25)
print(p)
}
# gse39582_Chen_MS4A4A.pdf
SurvivalByGene(gse_risk_norm,"MS4A4A")
# gse39582_Chen_C1orf162.pdf
SurvivalByGene(gse_risk_norm,"C1orf162")
# gse39582_Huang_MPP2.pdf
SurvivalByGene(gse_risk_norm,"MPP2")
# gse39582_Huang_PLEKHA8P1.pdf
SurvivalByGene(gse_risk_norm,"PLEKHA8P1")
# gse39582_Sun_RiskScore.pdf
## SUN
gse_risk_norm$RiskScoreSun<-0.1707*gse_risk_norm[,"FBXO17"]+0.1515*gse_risk_norm[,"PPARGC1A"]
SurvivalByGene(gse_risk_norm,"RiskScoreSun")
rm(exp)
rm(gse39582)
save.image("gse39582_CompareWithOthers.RData")
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/Rhistory/validation_gse39582_add.Rhistory")
