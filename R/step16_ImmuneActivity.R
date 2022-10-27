setwd("~/Desktop/Pro_Mine/Paper_TCGA_2/data_in")
library(utils)
library(estimate)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)    #need R 3.6
library(pheatmap)
library(ggpubr)
===================
# load gene sets for immune activity
sel_gmt=read.table("immune_features.txt",header = T,na.strings = c(" "),sep = "\t",stringsAsFactors = F)
sets=as.list(sel_gmt)
sets=lapply(sets, function(x) x[!is.na(x)])
sets<-lapply(sets, function(x) x[!x %in% ""])
# load expression FPKM
inputFile="symbol.txt"
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
# run GSVA
gsva_matrix<- gsva(as.matrix(data), sets,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
gsva_matrix1<- t(scale(t(gsva_matrix)))
normalization<-function(x){
return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1)
score.gsva.metab=as.data.frame(t(nor_gsva_matrix1))
# immune_activity_allSamples.pdf
nor_gsva_matrix1_tumor<-nor_gsva_matrix1[,40:ncol(nor_gsva_matrix1)]
p<-pheatmap(nor_gsva_matrix1_tumor,scale="row",show_colnames=F, show_rownames=T, cluster_cols=T, cluster_rows=F,cex=1, clustering_distance_rows="euclidean", cex=1,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,cutree_col = 2)
p
write.table(nor_gsva_matrix1_tumor,file="nor_gsva_matrix1_tumor.txt",sep="\t",row.names=F,quote=F)
# vln plot
cluster = cutree(p$tree_col,k=2)
table(cluster)
nor_gsva_matrix1_tumorVln<-colSums(nor_gsva_matrix1_tumor)
nor_gsva_matrix1_tumorVln<-as.data.frame(nor_gsva_matrix1_tumorVln)
cluster<-as.data.frame(cluster)
nor_gsva_matrix1_tumorVln<-cbind(nor_gsva_matrix1_tumorVln,cluster)
# immune_activity_vln_allSamples.pdf
ggplot(nor_gsva_matrix1_tumorVln, aes(factor(cluster), nor_gsva_matrix1_tumorVln))+ geom_violin(aes(fill = factor(cluster)))+xlab("cluster")+ylab("sum of immune enrichment score")

nor_gsva_matrix1_tumor_sum<-as.data.frame(t(nor_gsva_matrix1_tumor))
table(row.names(cluster)==row.names(nor_gsva_matrix1_tumor_sum))
TRUE 
 398 
nor_gsva_matrix1_tumor_sum<-cbind(nor_gsva_matrix1_tumor_sum,nor_gsva_matrix1_tumorVln)
write.table(nor_gsva_matrix1_tumor_sum,file="nor_gsva_matrix1_tumor_sum.txt",sep="\t",row.names=F,quote=F)
risk=read.table("risk.txt",sep="\t",header=T,row.names=1,check.names=F)
nor_gsva_matrix1_tumor_sum$SampleID=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",row.names(nor_gsva_matrix1_tumor_sum))
nor_gsva_matrix1_tumor_sum_risk334<-nor_gsva_matrix1_tumor_sum[match(row.names(risk),nor_gsva_matrix1_tumor_sum$SampleID),]
all.equal(row.names(risk),nor_gsva_matrix1_tumor_sum_risk334$SampleID)
[1] TRUE
nor_gsva_matrix1_tumor_sum_risk334<-cbind(nor_gsva_matrix1_tumor_sum_risk334,risk)
library(ggpubr)
nor_gsva_matrix1_tumor_sum_risk334$risk<-factor(nor_gsva_matrix1_tumor_sum_risk334$risk,levels = c("low","high"))
# immune_activity.pdf
pdf("immune_activity.pdf")
for (i in 1:30) {
p<-ggboxplot(nor_gsva_matrix1_tumor_sum_risk334, x = "risk", y = colnames(nor_gsva_matrix1_tumor_sum_risk334)[i],
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab(colnames(nor_gsva_matrix1_tumor_sum_risk334)[i])
print(p)
}
dev.off()

nor_gsva_matrix1_tumor_sum_risk334$riskScore_log2<-log2(nor_gsva_matrix1_tumor_sum_risk334$riskScore)


# immune_activity_riskscorelog2.pdf
ggboxplot(nor_gsva_matrix1_tumor_sum_risk334, x = "cluster", y = "riskScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("Risk score")
write.table(nor_gsva_matrix1_tumor_sum_risk334,file="nor_gsva_matrix1_tumor_sum_risk334.txt",sep="\t",row.names=F,quote=F)

nor_gsva_matrix1_tumor_sum_risk332<-nor_gsva_matrix1_tumor_sum_risk334[nor_gsva_matrix1_tumor_sum_risk334$riskScore<30,]

#### Boxplot
# immune_activity_riskscore_all334_yintercept4.pdf
ggboxplot(nor_gsva_matrix1_tumor_sum_risk334, x = "cluster", y = "riskScore",
color = "cluster", palette = "jco",shape=NA,
add = "jitter")+ stat_compare_means(method = "wilcox.test")+xlab("cluster")+ylab("Risk score")+geom_hline(yintercept=4)
dim(nor_gsva_matrix1_tumor_sum_risk334[nor_gsva_matrix1_tumor_sum_risk334$riskScore<4,])
nor_gsva_matrix1_tumor_sum_risk290<-nor_gsva_matrix1_tumor_sum_risk334[nor_gsva_matrix1_tumor_sum_risk334$riskScore<4,]
# immune_activity_riskscore_290.pdf
ggboxplot(nor_gsva_matrix1_tumor_sum_risk290, x = "cluster", y = "riskScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "wilcox.test")+xlab("cluster")+ylab("Risk score")
# immune_activity_riskscorelog2_all334.pdf
ggboxplot(nor_gsva_matrix1_tumor_sum_risk334, x = "cluster", y = "riskScore_log2",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "wilcox.test")+xlab("cluster")+ylab("log2(Risk score)")

#### Scatterplot
### Immune_activity_scatterplot_log2score_334.pdf
ggscatter(nor_gsva_matrix1_tumor_sum_risk334, x = "riskScore_log2", y = "nor_gsva_matrix1_tumorVln",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "log2 risk score", ylab = "sum of immune enrichment score")
### Immune_activity_scatterplot_score_334.pdf
ggscatter(nor_gsva_matrix1_tumor_sum_risk334, x = "riskScore", y = "nor_gsva_matrix1_tumorVln",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "risk score", ylab = "sum of immune enrichment score")
### Immune_activity_scatterplot_log2score_290.pdf
ggscatter(nor_gsva_matrix1_tumor_sum_risk290, x = "riskScore_log2", y = "nor_gsva_matrix1_tumorVln",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "log2 risk score", ylab = "sum of immune enrichment score")
### Immune_activity_scatterplot_score_290.pdf
ggscatter(nor_gsva_matrix1_tumor_sum_risk290, x = "riskScore", y = "nor_gsva_matrix1_tumorVln",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "risk score", ylab = "sum of immune enrichment score")



#### Dotheatmap
# Immune_activity_correlation_circle.pdf
library(corrplot)
## log2 risk score
M<-cor(nor_gsva_matrix1_tumor_sum_risk334[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk334))])
corrplot(M, method="circle")
M<-cor(nor_gsva_matrix1_tumor_sum_risk290[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk290))])
corrplot(M, method="circle")

#### Risk score
## Immune_activity_correlation_circle_290.pdf
M<-cor(nor_gsva_matrix1_tumor_sum_risk290[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk290)-2)])
corrplot(M, method="circle")
## Immune_activity_correlation_circle_334.pdf
M<-cor(nor_gsva_matrix1_tumor_sum_risk334[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk334)-2)])
corrplot(M, method="circle")

#### Dotheatmap upper 
#### Immune_activity_correlation_circle_2.pdf
cor.mtest <- function(mat, ...) {
mat <- as.matrix(mat)
n <- ncol(mat)
p.mat<- matrix(NA, n, n)
diag(p.mat) <- 0
for (i in 1:(n - 1)) {
for (j in (i + 1):n) {
tmp <- cor.test(mat[, i], mat[, j], ...)
p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
}
}
colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
p.mat
}


## Immune_activity_correlation_circle_2_290.pdf
p.mat <- cor.mtest(nor_gsva_matrix1_tumor_sum_risk290[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk290)-2)])
M<-cor(nor_gsva_matrix1_tumor_sum_risk290[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk290)-2)])
corrplot(M, type="upper", order="hclust", p.mat = p.mat, sig.level = 0.01, insig = "blank")
## Immune_activity_correlation_circle_2_334.pdf
p.mat <- cor.mtest(nor_gsva_matrix1_tumor_sum_risk334[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk334)-2)])
M<-cor(nor_gsva_matrix1_tumor_sum_risk334[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk334)-2)])
corrplot(M, type="upper", order="hclust", p.mat = p.mat, sig.level = 0.01, insig = "blank")



# Immune_activity_correlation_circle_4_290.pdf
M<-cor(nor_gsva_matrix1_tumor_sum_risk290[,c(1:29,ncol(nor_gsva_matrix1_tumor_sum_risk290))])
p.mat <- cor.mtest(nor_gsva_matrix1_tumor_sum_risk290[,c(1:29,ncol(nor_gsva_matrix1_tumor_sum_risk290)-2)])
corrplot(M, type="upper", p.mat = p.mat, sig.level = 0.01, insig = "blank")
# Immune_activity_correlation_circle_4_334.pdf
M<-cor(nor_gsva_matrix1_tumor_sum_risk334[,c(1:29,ncol(nor_gsva_matrix1_tumor_sum_risk334))])
p.mat <- cor.mtest(nor_gsva_matrix1_tumor_sum_risk334[,c(1:29,ncol(nor_gsva_matrix1_tumor_sum_risk334)-2)])
corrplot(M, type="upper", p.mat = p.mat, sig.level = 0.01, insig = "blank")

# heatmpa_imune_activity_RiskGenesExp_334.pdf
test<-nor_gsva_matrix1_tumor_sum_risk334[order(nor_gsva_matrix1_tumor_sum_risk334$cluster),]
test2<-t(test[,35:52])
anno_col<-as.data.frame(test$cluster)
row.names(anno_col)<-row.names(test)
colnames(anno_col)<-"cluster"
anno_col$cluster<-factor(anno_col$cluster,levels = c("1","2"))
pheatmap(test2,scale="row",cluster_cols = F,annotation_col = anno_col,show_colnames = F,color = colorRampPalette(c("blue", "white", "red"))(50))
write.table(test2,file="immune_activity_heatmap_test2_334.txt",sep="\t",row.names=F,quote=F)

# heatmpa_imune_activity_RiskGenesExp_290.pdf
test<-nor_gsva_matrix1_tumor_sum_risk290[order(nor_gsva_matrix1_tumor_sum_risk290$cluster),]
test2<-t(test[,35:52])
anno_col<-as.data.frame(test$cluster)
row.names(anno_col)<-row.names(test)
colnames(anno_col)<-"cluster"
anno_col$cluster<-factor(anno_col$cluster,levels = c("1","2"))
pheatmap(test2,scale="row",cluster_cols = F,annotation_col = anno_col,show_colnames = F,color = colorRampPalette(c("blue", "white", "red"))(50))
write.table(test2,file="immune_activity_heatmap_test2_290.txt",sep="\t",row.names=F,quote=F)


===================
##### cibersortx
cibersortx<-read.csv("CIBERSORTx_Job58_Results.csv",head=T,stringsAsFactors = F)
cibersortx$SampleID=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",cibersortx$Mixture)
cibersortx_risk<-cibersortx[match(row.names(risk),cibersortx$SampleID),]
cibersortx_risk<-cbind(cibersortx_risk,risk)
cibersortx_risk$riskScore_log2<-log2(cibersortx_risk$riskScore)
# correaltion_immunecells_riskscorelog2.pdf
pdf("correaltion_immunecells_riskscorelog2.pdf")
for (i in 2:23) {
p<-ggscatter(cibersortx_risk, x = "riskScore_log2", y = colnames(cibersortx_risk)[i],add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",xlab = "log2 risk score", ylab = colnames(cibersortx_risk)[i])
print(p)
}
dev.off()

# correaltion_immunecells_riskscorelog2_290.pdf
pdf("correaltion_immunecells_riskscorelog2_290.pdf")
for (i in 2:23) {
p<-ggscatter(cibersortx_risk_290, x = "riskScore_log2", y = colnames(cibersortx_risk)[i],add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",xlab = "log2 risk score", ylab = colnames(cibersortx_risk)[i])
print(p)
}
dev.off()

# correaltion_immunecells_riskscore_290.pdf
pdf("correaltion_immunecells_riskscore_290.pdf")
for (i in 2:23) {
p<-ggscatter(cibersortx_risk_290, x = "riskScore", y = colnames(cibersortx_risk)[i],add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",xlab = "risk score", ylab = colnames(cibersortx_risk)[i])
print(p)
}
dev.off()


#### Cell types from cibersortx
## 334 samples
cibersortx_risk<-cibersortx_risk[order(cibersortx_risk$risk),]
cibersortx_risk$risk<-factor(cibersortx_risk$risk,levels = c("low","high"))
cibersortx_risk<-cibersortx_risk[order(cibersortx_risk$risk),]
an_co<-data.frame("risk"=cibersortx_risk$risk)
row.names(an_co)<-cibersortx_risk$Mixture
cibersortx_risk_heatmap<-cibersortx_risk[,1:23]
row.names(cibersortx_risk_heatmap)<-cibersortx_risk_heatmap$Mixture
cibersortx_risk_heatmap<-cibersortx_risk_heatmap[,-1]
cibersortx_risk_heatmap<-t(cibersortx_risk_heatmap)
# Heatmap_riskgroup_immunecell.pdf
pheatmap(cibersortx_risk_heatmap,scale="row",cluster_cols = F,annotation_col = an_co,show_colnames = F,color = colorRampPalette(c("blue", "white", "red"))(50))
# boxplot_riskgroup_immunecell.pdf
pdf("boxplot_riskgroup_immunecell.pdf")
for (i in 2:23) {
p<-ggboxplot(cibersortx_risk, x = "risk", y = colnames(cibersortx_risk)[i],
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab(colnames(cibersortx_risk)[i])
print(p)
}
dev.off()
write.table(cibersortx_risk,file="cibersortx_risk.txt",sep="\t",row.names=F,quote=F)

## 290 samples
cibersortx_risk_290<-cibersortx_risk_290[order(cibersortx_risk_290$risk),]
cibersortx_risk_290$risk<-factor(cibersortx_risk_290$risk,levels = c("low","high"))
cibersortx_risk_290<-cibersortx_risk_290[order(cibersortx_risk_290$risk),]
an_co<-data.frame("risk"=cibersortx_risk_290$risk)
row.names(an_co)<-cibersortx_risk_290$Mixture
cibersortx_risk_290_heatmap<-cibersortx_risk_290[,1:23]
row.names(cibersortx_risk_290_heatmap)<-cibersortx_risk_290_heatmap$Mixture
cibersortx_risk_290_heatmap<-cibersortx_risk_290_heatmap[,-1]
cibersortx_risk_290_heatmap<-t(cibersortx_risk_290_heatmap)
# Heatmap_riskgroup_immunecell_290.pdf
pheatmap(cibersortx_risk_290_heatmap,scale="row",cluster_cols = F,annotation_col = an_co,show_colnames = F,color = colorRampPalette(c("blue", "white", "red"))(50))
# boxplot_riskgroup_immunecell_290.pdf
pdf("boxplot_riskgroup_immunecell_290.pdf")
for (i in 2:23) {
p<-ggboxplot(cibersortx_risk_290, x = "risk", y = colnames(cibersortx_risk_290)[i],
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab(colnames(cibersortx_risk_290)[i])
print(p)
}
dev.off()
write.table(cibersortx_risk_290,file="cibersortx_risk_290.txt",sep="\t",row.names=F,quote=F)

===================
##### ESTIMATE score
estimate=read.table("COAD_estimate_score.gct.txt",sep="\t",header=T,check.names=F,stringsAsFactors = F)
estimate<-estimate[,42:ncol(estimate)]
colnames(estimate)<-gsub("\\.","-",colnames(estimate))
colnames(estimate)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(estimate))
estimate<-estimate[,row.names(risk)]
row.names(estimate)<-c("StromalScore","ImmuneScore","ESTIMATEScore")
estimate<-t(estimate)
all.equal(row.names(estimate),nor_gsva_matrix1_tumor_sum_risk334$SampleID)
activity_risk_estimate_334<-cbind(nor_gsva_matrix1_tumor_sum_risk334,estimate)
library(ggpubr)
pdf("estimate_334.pdf")
ggboxplot(activity_risk_estimate_334, x = "risk", y = "ESTIMATEScore",
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab("ESTIMATEScore")
ggboxplot(activity_risk_estimate_334, x = "risk", y = "StromalScore",
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab("StromalScore")
ggboxplot(activity_risk_estimate_334, x = "risk", y = "ImmuneScore",
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab("ImmuneScore")
ggboxplot(activity_risk_estimate_334, x = "cluster", y = "ImmuneScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("ImmuneScore")
ggboxplot(activity_risk_estimate_334, x = "cluster", y = "StromalScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("StromalScore")
ggboxplot(activity_risk_estimate_334, x = "cluster", y = "ESTIMATEScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("ESTIMATEScore")
activity_risk_estimate_334$TumorScore<-max(activity_risk_estimate_334$ESTIMATEScore)-activity_risk_estimate_334$ESTIMATEScore
ggboxplot(activity_risk_estimate_334, x = "cluster", y = "TumorScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("TumorScore")
ggboxplot(activity_risk_estimate_334, x = "risk", y = "TumorScore",
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab("TumorScore")
ggscatter(activity_risk_estimate_334, x = "riskScore_log2", y = "TumorScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "log2 risk score", ylab = "TumorScore")
ggscatter(activity_risk_estimate_334, x = "riskScore_log2", y = "ESTIMATEScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "log2 risk score", ylab = "ESTIMATEScore")
ggscatter(activity_risk_estimate_334, x = "riskScore_log2", y = "ImmuneScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "log2 risk score", ylab = "ImmuneScore")
ggscatter(activity_risk_estimate_334, x = "nor_gsva_matrix1_tumorVln", y = "ImmuneScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "sum of immune enrichment score", ylab = "ImmuneScore")
ggscatter(activity_risk_estimate_334, x = "nor_gsva_matrix1_tumorVln", y = "ESTIMATEScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "sum of immune enrichment score", ylab = "ESTIMATEScore")
ggscatter(activity_risk_estimate_334, x = "nor_gsva_matrix1_tumorVln", y = "TumorScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "sum of immune enrichment score", ylab = "TumorScore")
ggscatter(activity_risk_estimate_334, x = "nor_gsva_matrix1_tumorVln", y = "StromalScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "sum of immune enrichment score", ylab = "StromalScore")
dev.off()
write.table(activity_risk_estimate_334,file="activity_risk_estimate_334.txt",sep="\t",row.names=F,quote=F)
activity_risk_estimate_290<-activity_risk_estimate_334[activity_risk_estimate_334$riskScore<4,]
dim(activity_risk_estimate_290)

pdf("estimate_290.pdf")
ggboxplot(activity_risk_estimate_290, x = "risk", y = "ESTIMATEScore",
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab("ESTIMATEScore")
ggboxplot(activity_risk_estimate_290, x = "risk", y = "StromalScore",
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab("StromalScore")
ggboxplot(activity_risk_estimate_290, x = "risk", y = "ImmuneScore",
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab("ImmuneScore")
ggboxplot(activity_risk_estimate_290, x = "cluster", y = "ImmuneScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("ImmuneScore")
ggboxplot(activity_risk_estimate_290, x = "cluster", y = "StromalScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("StromalScore")
ggboxplot(activity_risk_estimate_290, x = "cluster", y = "ESTIMATEScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("ESTIMATEScore")
activity_risk_estimate_290$TumorScore<-max(activity_risk_estimate_290$ESTIMATEScore)-activity_risk_estimate_290$ESTIMATEScore
ggboxplot(activity_risk_estimate_290, x = "cluster", y = "TumorScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("TumorScore")
ggboxplot(activity_risk_estimate_290, x = "risk", y = "TumorScore",
color = "risk", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab("TumorScore")
ggscatter(activity_risk_estimate_290, x = "riskScore_log2", y = "TumorScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "log2 risk score", ylab = "TumorScore")
ggscatter(activity_risk_estimate_290, x = "riskScore_log2", y = "ESTIMATEScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "log2 risk score", ylab = "ESTIMATEScore")
ggscatter(activity_risk_estimate_290, x = "riskScore_log2", y = "ImmuneScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "log2 risk score", ylab = "ImmuneScore")
ggscatter(activity_risk_estimate_290, x = "nor_gsva_matrix1_tumorVln", y = "ImmuneScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "sum of immune enrichment score", ylab = "ImmuneScore")
ggscatter(activity_risk_estimate_290, x = "nor_gsva_matrix1_tumorVln", y = "ESTIMATEScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "sum of immune enrichment score", ylab = "ESTIMATEScore")
ggscatter(activity_risk_estimate_290, x = "nor_gsva_matrix1_tumorVln", y = "TumorScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "sum of immune enrichment score", ylab = "TumorScore")
ggscatter(activity_risk_estimate_290, x = "nor_gsva_matrix1_tumorVln", y = "StromalScore",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "sum of immune enrichment score", ylab = "StromalScore")
dev.off()
write.table(activity_risk_estimate_290,file="activity_risk_estimate_290.txt",sep="\t",row.names=F,quote=F)

# immune_activity_2_334.pdf
long<-reshape2::melt(nor_gsva_matrix1_tumor_sum_risk334[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk334)-1)],
id.vars = c('risk'),
variable.name='geneSets',
value.name='ssGSEAscores')
ggboxplot(long, x = "risk", y = "ssGSEAscores",
color = "risk", palette = "jco",
facet.by = "geneSets", short.panel.labs = FALSE,scales='free') + stat_compare_means(label = "p.format")
# immune_activity_2_290.pdf
long_290<-reshape2::melt(nor_gsva_matrix1_tumor_sum_risk290[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk290)-1)],
id.vars = c('risk'),
variable.name='geneSets',
value.name='ssGSEAscores')
ggboxplot(long_290, x = "risk", y = "ssGSEAscores",
color = "risk", palette = "jco",
facet.by = "geneSets", short.panel.labs = FALSE,scales='free') + stat_compare_means(label = "p.format")

# Clusters_riskGenes_334.pdf
clusters<-activity_risk_estimate_334[,c(31,35:52,54,56,57,59)]
clusters<-clusters[order(clusters$cluster),]
anno_clusters<-clusters[,c(1,20,21,22,23)]
anno_clusters<-anno_clusters[,c(5,4,3,2,1)]
anno_clusters$cluster<-factor(anno_clusters$cluster,levels = c("1","2"))
pheatmap(t(clusters[,c(2:19)]),scale="row",show_colnames=F, show_rownames=T, cluster_cols=F, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", cex=1,fontsize = 8,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,annotation_col = anno_clusters,color = colorRampPalette(c("blue", "white", "red"))(50))
# Clusters_ISR_Exp_MeanZscore_334.pdf
clusters_boxplot<-clusters[,c(1:19)]
clusters_boxplot$MeanZscore<-rowMeans(scale(clusters_boxplot[,2:ncol(clusters_boxplot)]))
clusters_boxplot$cluster<-factor(clusters_boxplot$cluster,levels = c("1","2"))
ggboxplot(clusters_boxplot, x = "cluster", y = "MeanZscore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")

# Clusters_riskGenes_290.pdf
clusters<-activity_risk_estimate_290[,c(31,35:52,54,56,57,59)]
clusters<-clusters[order(clusters$cluster),]
anno_clusters<-clusters[,c(1,20,21,22,23)]
anno_clusters<-anno_clusters[,c(5,4,3,2,1)]
anno_clusters$cluster<-factor(anno_clusters$cluster,levels = c("1","2"))
pheatmap(t(clusters[,c(2:19)]),scale="row",show_colnames=F, show_rownames=T, cluster_cols=F, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", cex=1,fontsize = 8,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,annotation_col = anno_clusters,color = colorRampPalette(c("blue", "white", "red"))(50))
# Clusters_ISR_Exp_MeanZscore_290.pdf
clusters_boxplot<-clusters[,c(1:19)]
clusters_boxplot$MeanZscore<-rowMeans(scale(clusters_boxplot[,2:ncol(clusters_boxplot)]))
clusters_boxplot$cluster<-factor(clusters_boxplot$cluster,levels = c("1","2"))
ggboxplot(clusters_boxplot, x = "cluster", y = "MeanZscore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")

===========
rm(rt)
rm(exp)
rm(data)
save.image("immuneActivity.RData")
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/Rhistory/pip2_ImmuneActivity.Rhistory")





# Clusters_riskGenes_290_redo.pdf  #by risk category
clusters<-activity_risk_estimate_290
clusters<-clusters[,c(31,35:57,59,61:64)]
clusters<-clusters[order(clusters$cluster),]
clusters[,c(1,24,25,26,27)]
anno_clusters<-clusters[,c(1,25,26,27,28,29)]
anno_clusters[,c(6,5,4,3,2,1)]
colnames(anno_clusters[,c(6,5,4,3,2,1)])
anno_clusters<-anno_clusters[,c(6,5,4,3,2,1)]
anno_clusters$cluster<-factor(anno_clusters$cluster,levels = c("1","2"))
pheatmap(t(clusters[,c(2:24)]),scale="row",show_colnames=F, show_rownames=T, cluster_cols=F, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", cex=1,fontsize = 8,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,annotation_col = anno_clusters,color = colorRampPalette(c("blue", "white", "red"))(50))

# Clusters_riskGenes_290_redo2.pdf   #by risk score
clusters<-activity_risk_estimate_290
clusters<-clusters[,c(31,35:57,58,61:64)]
clusters<-clusters[order(clusters$cluster),]
anno_clusters<-clusters[,c(1,25,26,27,28,29)]
colnames(anno_clusters[,c(6,5,4,3,2,1)])
anno_clusters<-anno_clusters[,c(6,5,4,3,2,1)]
anno_clusters$cluster<-factor(anno_clusters$cluster,levels = c("1","2"))
pheatmap(t(clusters[,c(2:24)]),scale="row",show_colnames=F, show_rownames=T, cluster_cols=F, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", cex=1,fontsize = 8,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,annotation_col = anno_clusters,color = colorRampPalette(c("blue", "white", "red"))(50))
deg_cluster<-clusters %>% dplyr::select(-cluster,cluster)
class(deg_cluster$riskScore)

## outTab_2 loose
deg_cluster<-clusters %>% dplyr::select(-cluster,cluster)
outTab_2=data.frame()
for(i in colnames(deg_cluster[,1:ncol(deg_cluster)-1])){
geneName=i
rt=rbind(expression=deg_cluster[,i],grade=deg_cluster[,ncol(deg_cluster)])
rt=as.matrix(t(rt))
wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
conGeneMeans=mean(deg_cluster[deg_cluster$cluster==1,i])
treatGeneMeans=mean(deg_cluster[deg_cluster$cluster==2,i])
logFC=log2(treatGeneMeans)-log2(conGeneMeans)
pvalue=wilcoxTest$p.value
conMed=median(deg_cluster[deg_cluster$cluster==1,i])
treatMed=median(deg_cluster[deg_cluster$cluster==2,i])
diffMed=treatMed-conMed
outTab_2=rbind(outTab_2,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
}
pValue=outTab_2[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab_2=cbind(outTab_2,fdr=fdr)
outTab_2[outTab_2$fdr<0.05,]
write.csv(outTab_2,"./DEGs_ByCluster_2.csv")
outTab_2[outTab_2$fdr>=0.05,]
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/Rhistory/pip3_ImmuneActivity.Rhistory")


##### Define outliers with grubbs.test

#install.packages("outliers")
library(outliers)
outliers.test.334<-nor_gsva_matrix1_tumor_sum_risk334$riskScore

# check minium if it is outlier
outliers.test.334<-nor_gsva_matrix1_tumor_sum_risk334$riskScore
for (i in 1:length(outliers.test.334)) {
gtest<-grubbs.test(outliers.test.334,opposite=TRUE)
if (gtest$p.value<0.05) {
print(gtest$alternative)
print(min(outliers.test.334))
outliers.test.334<-outliers.test.334[outliers.test.334!=min(outliers.test.334)]
} else {
break
}
}
# no outliers in minum values

# check maximum if it is outlier
outliers.test.334<-nor_gsva_matrix1_tumor_sum_risk334$riskScore
for (i in 1:length(outliers.test.334)) {
gtest<-grubbs.test(outliers.test.334)
if (gtest$p.value<0.05) {
print(gtest$alternative)
print(max(outliers.test.334))
outliers.test.334<-outliers.test.334[outliers.test.334!=max(outliers.test.334)]
} else {
break
}
}

outliers.test.313<-outliers.test.334
nor_gsva_matrix1_tumor_sum_risk313<-nor_gsva_matrix1_tumor_sum_risk334[nor_gsva_matrix1_tumor_sum_risk334$riskScore %in% outliers.test.313,]

####

# immune_activity_riskscore_313.pdf # this time change to t.text now because of normal distribtuion
ggboxplot(nor_gsva_matrix1_tumor_sum_risk313, x = "cluster", y = "riskScore",
color = "cluster", palette = "jco",
add = "jitter")+ stat_compare_means(method = "t.test")+xlab("cluster")+ylab("Risk score")

#activity_risk_estimate_334<-cbind(nor_gsva_matrix1_tumor_sum_risk334,estimate)
activity_risk_estimate_313<-activity_risk_estimate_334[activity_risk_estimate_334$riskScore<=max(nor_gsva_matrix1_tumor_sum_risk313$riskScore),]
# Clusters_riskGenes_313_redo.pdf  #by risk category 6.73*4.56
clusters<-activity_risk_estimate_313
clusters<-clusters[,c(31,35:57,59,61:64)]
clusters<-clusters[order(clusters$cluster),]
clusters[,c(1,24,25,26,27)]
anno_clusters<-clusters[,c(1,25,26,27,28,29)]
anno_clusters[,c(6,5,4,3,2,1)]
colnames(anno_clusters[,c(6,5,4,3,2,1)])
anno_clusters<-anno_clusters[,c(6,5,4,3,2,1)]
anno_clusters$cluster<-factor(anno_clusters$cluster,levels = c("1","2"))
pheatmap(t(clusters[,c(2:24)]),scale="row",show_colnames=F, show_rownames=T, cluster_cols=F, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", cex=1,fontsize = 8,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,annotation_col = anno_clusters,color = colorRampPalette(c("blue", "white", "red"))(50))
# Clusters_riskGenes_313_redo2.pdf   #by risk score  6.73*4.56
clusters<-activity_risk_estimate_313
clusters<-clusters[,c(31,35:57,58,61:64)]
clusters<-clusters[order(clusters$cluster),]
anno_clusters<-clusters[,c(1,25,26,27,28,29)]
colnames(anno_clusters[,c(6,5,4,3,2,1)])
anno_clusters<-anno_clusters[,c(6,5,4,3,2,1)]
anno_clusters$cluster<-factor(anno_clusters$cluster,levels = c("1","2"))
pheatmap(t(clusters[,c(2:24)]),scale="row",show_colnames=F, show_rownames=T, cluster_cols=F, cluster_rows=T,cex=1, clustering_distance_rows="euclidean", cex=1,fontsize = 8,clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE,annotation_col = anno_clusters,color = colorRampPalette(c("blue", "white", "red"))(50))
deg_cluster<-clusters %>% dplyr::select(-cluster,cluster)
class(deg_cluster$riskScore)

### Immune_activity_scatterplot_score_313.pdf 4.34*4.56
ggscatter(nor_gsva_matrix1_tumor_sum_risk313, x = "riskScore", y = "nor_gsva_matrix1_tumorVln",
add = "reg.line", conf.int = TRUE,
cor.coef = TRUE, cor.method = "pearson",
xlab = "risk score", ylab = "sum of immune enrichment score")
## Immune_activity_correlation_circle_2_313.pdf 12*12
p.mat <- cor.mtest(nor_gsva_matrix1_tumor_sum_risk313[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk313)-2)])
M<-cor(nor_gsva_matrix1_tumor_sum_risk313[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk313)-2)])
corrplot(M, type="upper", order="hclust", p.mat = p.mat, sig.level = 0.01, insig = "blank")
## outTab_2 loose
deg_cluster<-clusters %>% dplyr::select(-cluster,cluster)
outTab_2=data.frame()
for(i in colnames(deg_cluster[,1:ncol(deg_cluster)-1])){
geneName=i
rt=rbind(expression=deg_cluster[,i],grade=deg_cluster[,ncol(deg_cluster)])
rt=as.matrix(t(rt))
wilcoxTest<-wilcox.test(expression ~ grade, data=rt)
conGeneMeans=mean(deg_cluster[deg_cluster$cluster==1,i])
treatGeneMeans=mean(deg_cluster[deg_cluster$cluster==2,i])
logFC=log2(treatGeneMeans)-log2(conGeneMeans)
pvalue=wilcoxTest$p.value
conMed=median(deg_cluster[deg_cluster$cluster==1,i])
treatMed=median(deg_cluster[deg_cluster$cluster==2,i])
diffMed=treatMed-conMed
outTab_2=rbind(outTab_2,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
}
pValue=outTab_2[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)),method="fdr")
outTab_2=cbind(outTab_2,fdr=fdr)
outTab_2[outTab_2$fdr<0.05,]
write.csv(outTab_2,"./DEGs_ByCluster_2_313samples.csv")
outTab_2[outTab_2$fdr>=0.05,]
save.image("immuneActivity.RData")
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/Rhistory/pip4_ImmuneActivity.Rhistory")
# immune_activity_2_313.pdf # 20*20 ## "nor_gsva_matrix1_tumorVln" is SIES
long_313<-reshape2::melt(nor_gsva_matrix1_tumor_sum_risk313[,c(1:30,ncol(nor_gsva_matrix1_tumor_sum_risk313)-1)],
id.vars = c('risk'),
variable.name='geneSets',
value.name='ssGSEAscores')
ggboxplot(long_313, x = "risk", y = "ssGSEAscores",
color = "risk", palette = "jco",
facet.by = "geneSets", short.panel.labs = FALSE,scales='free') + stat_compare_means(label = "p.format")
# boxplot_riskgroup_immunecell_313.pdf
pdf("boxplot_riskgroup_immunecell_313.pdf")
for (i in 2:23) {
    p<-ggboxplot(cibersortx_risk_313, x = "risk", y = colnames(cibersortx_risk_313)[i],
                 color = "risk", palette = "jco",
                 add = "jitter")+ stat_compare_means(method = "t.test")+xlab("risk")+ylab(colnames(cibersortx_risk_313)[i])
    print(p)
}
dev.off()


# Immune_activity_boxplot_score_313.pdf ## "nor_gsva_matrix1_tumorVln" is SIES
ggboxplot(nor_gsva_matrix1_tumor_sum_risk313, x = "risk", y = "nor_gsva_matrix1_tumorVln",color="risk", palette="jco")+stat_compare_means(method = "t.test")+xlab("risk group")+ylab("sum of immune enrichment score")

# immune_activity_2_313_Ttest.pdf ## "nor_gsva_matrix1_tumorVln" is SIES
ggboxplot(long_313, x = "risk", y = "ssGSEAscores",
          color = "risk", palette = "jco",
          facet.by = "geneSets", short.panel.labs = FALSE,scales='free') + stat_compare_means(label = "p.format",method = "t.test")
