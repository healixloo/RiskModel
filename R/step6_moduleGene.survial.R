
#install.packages("survival")
#install.packages("survminer")

setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")      
library(survival)
library("survminer")
rt=read.table("risk.txt",header=T,sep="\t")
diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)

# survival curve
pdf(file="survival.pdf",onefile = FALSE,
       width = 5.5,             
       height =5)             
ggsurvplot(fit, 
           data=rt,
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
dev.off()
sink("fit.txt")
summary(fit)
sink()

