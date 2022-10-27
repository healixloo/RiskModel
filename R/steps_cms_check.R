library(ggpubr)
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")
rt=read.table("risk.txt",header=T,sep="\t")
clinical_cms<-read.table("./clinical_molecular_public_all.txt",header = T,stringsAsFactors = F, check.names = F)
annotation_cms<-read.table("cms_labels_public_all.txt",header = T,stringsAsFactors = F, check.names = F)
clinical_cbio<-read.csv("cbioportal_combined_study_clinical_data.tsv",header = T,stringsAsFactors = F, check.names = F,sep = '\t')
rt_annotation<-merge(rt,clinical_cms,by.x = "id",by.y="sample",all.x=T)
rt_annotation<-merge(rt_annotation,annotation_cms,by.x = "id",by.y="sample",all.x=T)
table(rt_annotation$cms_label==rt_annotation$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples)

colnames(clinical_cbio)<-paste("cBio_",colnames(clinical_cbio),sep="")
colnames(clinical_cbio)<-gsub(" ","_",colnames(clinical_cbio))
# rt_annotation<-merge(rt_annotation,clinical_cbio,by.x = "id",by.y="cBio_Patient_ID",all.x=T)
write.csv(rt_annotation,"rt_annotation_cms.csv")

rt_annotation$log2_riskScore<-log2(rt_annotation$riskScore)
## cms_vs_log2riskScore.pdf
rt_annotation$cms_label<-factor(rt_annotation$cms_label,levels = c("CMS1","CMS2","CMS3","CMS4","NOLBL"))
ggboxplot(data=subset(rt_annotation, !is.na(cms_label)), x = "cms_label", y = "log2_riskScore",size=0.5,
color = "cms_label", palette = "jco",
add = "jitter")+stat_compare_means(comparisons = list( c("CMS1", "CMS4"), c("CMS2", "CMS4"), c("CMS3", "CMS4"),c("CMS4","NOLBL") ),method = "t.test")+scale_colour_grey()
> table(rt_annotation$cms_label)
 CMS1  CMS2  CMS3  CMS4 NOLBL 
   40   121    37    75    32 
> table(rt_annotation$risk,rt_annotation$cms_label)
      
       CMS1 CMS2 CMS3 CMS4 NOLBL
  high   14   57   15   48    15
  low    26   64   22   27    17
## cms_vs_NCOA7.pdf
ggboxplot(data=subset(rt_annotation, !is.na(cms_label)), x = "cms_label", y = "NCOA7",size=0.5,
color = "cms_label", palette = "jco",
add = "jitter")+stat_compare_means(comparisons = list( c("CMS1", "CMS2"), c("CMS1", "CMS3"), c("CMS1", "CMS4"),c("CMS1","NOLBL") ),method = "t.test")+scale_colour_grey()

rt_annotation_df<-as.data.frame(table(rt_annotation$risk,rt_annotation$cms_label))
library(dplyr)
rt_annotation_df<-rt_annotation_df %>% group_by(Var1) %>% mutate(Ratio=Freq/sum(Freq))
rt_annotation_df$Var1<-factor(rt_annotation_df$Var1,levels = c("low","high"))
## simple check by stacked barplot
ggplot(rt_annotation_df, aes(fill=Var2, y=Freq, x=Var1)) +
geom_bar(position="fill", stat="identity")
## CMS_composistion.pdf
ggplot(rt_annotation_df, aes(fill=Var2, y=Ratio, x=Var1)) +
geom_bar(position="stack", stat="identity")+geom_text(aes(label = round(Ratio,digits=3)), position = position_stack(vjust = 0.5), size = 4)+scale_fill_manual(values = c("#fd9326","#1ca4fc","#ed62a7","#29af1d","#ababaf"))


#### msi_vs_log2riskScore.pdf
ggboxplot(data=subset(rt_annotation, !is.na(msi)), x = "msi", y = "log2_riskScore",size=0.5,
color = "msi", palette = "jco",
add = "jitter")+stat_compare_means(method = "t.test")+scale_colour_grey()
table(rt_annotation$msi)
msi mss 
 42 253 

### kras_vs_log2riskScore.pdf 
ggboxplot(data=subset(rt_annotation, !is.na(kras_mut)), x = "kras_mut", y = "log2_riskScore",size=0.5,
color = "kras_mut", palette = "jco",
add = "jitter")+stat_compare_means(method = "t.test")+scale_colour_grey()
table(rt_annotation$kras_mut)
  0   1 
126  77 
### braf_vs_log2riskScore.pdf
ggboxplot(data=subset(rt_annotation, !is.na(braf_mut)), x = "braf_mut", y = "log2_riskScore",size=0.5,
color = "braf_mut", palette = "jco",
add = "jitter")+stat_compare_means(method = "t.test")+scale_colour_grey()
table(rt_annotation$braf_mut)
  0   1 
184  19 
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/steps_cms_check.R")
### Kras_NCOA7.pdf
ggboxplot(data=subset(rt_annotation, !is.na(kras_mut)), x = "kras_mut", y = "NCOA7",size=0.5,
          color = "kras_mut", palette = "jco",
          add = "jitter")+stat_compare_means(method = "t.test")+scale_colour_grey()

### Braf_NCOA7.pdf
ggboxplot(data=subset(rt_annotation, !is.na(braf_mut)), x = "braf_mut", y = "NCOA7",size=0.5,
          color = "braf_mut", palette = "jco",
          add = "jitter")+stat_compare_means(method = "t.test")+scale_colour_grey()

          
## Mutation load
# MutationCount_vs_log2riskScore.pdf
rt_annotation<-merge(rt_annotation,clinical_cbio,by.x = "id",by.y="cBio_Patient_ID",all.x=T)
rt_annotation_atlas<-rt_annotation[rt_annotation$cBio_Study_ID=="coadread_tcga_pan_can_atlas_2018",]
ggscatter(rt_annotation_atlas, x = "cBio_Mutation_Count", y = "log2_riskScore",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) + stat_cor(method = "pearson")
ggscatter(rt_annotation_atlas, x = "cBio_Mutation_Count", y = "log2_riskScore",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) + stat_cor(method = "pearson")
nrow(rt_annotation_atlas)
 318
write.csv(rt_annotation,"rt_annotation_cBio.csv")
write.csv(rt_annotation_atlas,"rt_annotation_atlas.csv")

### check KRAS and BRAF expression vs riskScore
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")
library(limma)
inputFile="symbol.txt"
fdrFilter=0.05
logFCfilter=1
#conNum=39
#treatNum=398
# Read in data, Average gene when genes were duplicated in matrix
outTab=data.frame()
#grade=c(rep(1,conNum),rep(2,treatNum))
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data_tumor<-data[,c(40:ncol(data))]
risk=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)
grade<-as.character(risk[substring(colnames(data_tumor),1,12),"risk"])
data_tumor_targets<-data_tumor[c("KRAS","BRAF"),]
data_tumor_targets<-as.data.frame(t(data_tumor_targets))
data_tumor_targets$sID<-substring(row.names(data_tumor_targets),1,12)
data_tumor_targets<-cbind(data_tumor_targets,risk[data_tumor_targets$sID,])
data_tumor_targets$log2_riskScore<-log2(data_tumor_targets$riskScore)

library(ggpubr)
# risk_vs_expKRAS.pdf
ggboxplot(data=subset(data_tumor_targets, !is.na(risk)), x = "risk", y = "KRAS",size=0.5,
color = "risk", palette = "jco",
add = "jitter")+stat_compare_means(method = "t.test")+scale_colour_grey()
# risk_vs_expBRAF
ggboxplot(data=subset(data_tumor_targets, !is.na(risk)), x = "risk", y = "BRAF",size=0.5,
color = "risk", palette = "jco",
add = "jitter")+stat_compare_means(method = "t.test")+scale_colour_grey()
table(data_tumor_targets$risk)
high  low 
 178  172 
# expKRAS_vs_log2_riskScore.pdf expBRAF_vs_log2_riskScore.pdf
ggscatter(data_tumor_targets, x = "KRAS", y = "log2_riskScore",
add = "reg.line",  # Add regressin line
add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
conf.int = TRUE # Add confidence interval
) + stat_cor(method = "pearson")
ggscatter(data_tumor_targets, x = "BRAF", y = "log2_riskScore",
add = "reg.line",  # Add regressin line
add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
conf.int = TRUE # Add confidence interval
) + stat_cor(method = "pearson")


=====================================================================
library(ggpubr)
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")
rt=read.table("risk.txt",header=T,sep="\t")
#treatment_all<-read.table("./media-2.tsv",header = T,stringsAsFactors = F,check.names = F,sep = '\t')
#treatment_coad<-treatment_all[grepl("COAD",treatment_all$tumor),]
#treatment_coad$sid<-toupper(treatment_coad$patient_barcode)
#rtm<-merge(rt,treatment_coad,by.x = "id",by.y="sid",all.x = T)
#rtm<-rtm[!is.na(rtm$therapy_type),]
#rtm$log2_riskScore<-log2(rtm$riskScore)
#ggboxplot(rtm[rtm$tumor=="COAD",], x = "measure_of_response", y = "log2_riskScore",size=0.5,color = "measure_of_response", palette = "jco",add = "jitter")

treatment_filter<-read.table("./media-4.tsv",header = T,stringsAsFactors = F,check.names = F,sep = '\t')
rtm_filter<-merge(rt,treatment_filter,by.x = "id",by.y="patient_ID",all.x = T)
rtm_filter<-merge(rt,treatment_filter,by.x = "id",by.y="patient_ID",all.x = T)
rtm_filter_COAD<-rtm_filter[rtm_filter$tumor=="COAD" & (!is.na(rtm_filter$response)),]
rtm_filter_COAD$log2_riskScore<-log2(rtm_filter_COAD$riskScore)
rtm_filter_COAD$response<-factor(rtm_filter_COAD$response,levels = c("NON RESPONDER","RESPONDER"))

# Therapy_vs_riskScore.pdf
ggboxplot(rtm_filter_COAD, x = "response", y = "log2_riskScore",size=0.5,
color = "response", palette = "jco",
add = "jitter")+stat_compare_means(method = "t.test")+scale_color_manual(values = c("grey","grey0"))
> table(rtm_filter_COAD$response)
NON RESPONDER     RESPONDER 
            5            20 
> table(rtm_filter_COAD$treatment)
FLUOROURACIL+LEUCOVORIN+OXALIPLATIN 
                                 25 
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/steps_cms_check_pip4.Rhistory")

### NCOA7 Expression in tumor and normal
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")
rt=read.table("risk.txt",header=T,sep="\t")
tcga_exp=read.table("symbol.txt",sep="\t",header=T,check.names=F)
tcga_exp_t<-as.data.frame(t(tcga_exp))
colnames(tcga_exp_t)<-tcga_exp$ID
tcga_exp_t<-tcga_exp_t[-1,]
tcga_exp_t_NCOA7<-as.data.frame(tcga_exp_t[,"NCOA7"])
colnames(tcga_exp_t_NCOA7)<-"NCOA7"
tcga_exp_t_NCOA7$pid<-substring(row.names(tcga_exp_t_NCOA7),1,12)
tcga_exp_t_NCOA7$type<-substring(row.names(tcga_exp_t_NCOA7),14,15)
tcga_exp_t_NCOA7$types<-"tumor"
tcga_exp_t_NCOA7[tcga_exp_t_NCOA7$type=="11","types"]<-"normal"

rtx<-rt[,c("id","NCOA7","riskScore","risk")]
rty<-tcga_exp_t_NCOA7[tcga_exp_t_NCOA7$types=="normal",c("pid","NCOA7","types","types")]
colnames(rtx)<-c("pid","NCOA7","riskScore","types")
colnames(rty)[3]<-"riskScore"
colnames(rty)[4]<-"types"
rtz<-rbind(rtx,rty)
rtz$Type<-"tumor"
rtz[rtz$types=="normal","Type"]<-"normal"

rtz$NCOA7<-as.numeric(rtz$NCOA7)
rtz$NCOA7_log2<-log2(rtz$NCOA7)

rtz$types<-factor(rtz$types,levels = c("normal","low","high"))
rtz$Type<-factor(rtz$Type,levels = c("normal","tumor"))

## NCOA7_expLog2_All_1.pdf
ggboxplot(rtz, x = "types", y = "NCOA7_log2",size=0.5,
color = "types", palette = "jco",
add = "jitter")+stat_compare_means(comparisons = list(c("normal","low"),c("normal","high"),c("low","high")),method = "t.test")+scale_color_manual(values = c("grey","grey0","grey0"))
> table(rtz$types)
normal    low   high 
    39    167    167 
## NCOA7_expLog2_All_2.pdf
ggboxplot(rtz, x = "Type", y = "NCOA7_log2",size=0.5,
color = "Type", palette = "jco",
add = "jitter")+stat_compare_means(method = "t.test")+scale_color_manual(values = c("grey","grey0"))
> table(rtz$Type)
normal  tumor 
    39    334 

savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/steps_cms_check_pip5.Rhistory")




################# Additional check
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_in")
rt_annotation<-read.csv("rt_annotation_cms.csv",header = T,stringsAsFactors = F,check.names = F)
rt_annotation$log2_riskScore<-log2(rt_annotation$riskScore)
rt_annotation_2<-read.csv("rt_annotation_cBio.csv",header = T,stringsAsFactors = F,check.names = F)
library(ggpubr)
# cimp_vs_riskScore.pdf
ggboxplot(data=rt_annotation, x = "cimp", y = "log2_riskScore",size=0.5,
color = "cimp", palette = "jco",
add = "jitter")+stat_compare_means(comparisons = list( c("CIMP.Neg", "CIMP.Low"), c("CIMP.Neg", "CIMP.High"), c("CIMP.Low", "CIMP.High")),method = "t.test")+scale_colour_grey()
# copyNum_vs_riskScore.pdf
ggboxplot(data=rt_annotation_2, x = "cBio_Copy_Number", y = "log2_riskScore",size=0.5,
color = "cBio_Copy_Number", palette = "jco",
add = "jitter")+stat_compare_means(comparisons = list( c("0", "1")),method = "t.test")+scale_colour_grey()
# HyperMutated_vs_riskScore.pdf
rt_annotation_2$cBio_Hyper_mutated<-rt_annotation_2$`cBio_Hyper-mutated`
ggboxplot(data=rt_annotation_2, x = "cBio_Hyper_mutated", y = "log2_riskScore",size=0.5,
color = "cBio_Hyper_mutated", palette = "jco",
add = "jitter")+stat_compare_means(comparisons = list( c("0", "1")),method = "t.test")+scale_colour_grey()
# Methylation_vs_riskScore.pdf
ggboxplot(data=rt_annotation_2, x = "cBio_Methylation_Status", y = "log2_riskScore",size=0.5,
color = "cBio_Methylation_Status", palette = "jco",
add = "jitter")+stat_compare_means(comparisons = list( c("0", "1")),method = "t.test")+scale_colour_grey()
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/Rhistory/steps_cms_check_pip7.Rhistory")
