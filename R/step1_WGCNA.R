
setwd("~/Desktop/Pro_Mine/Paper_TCGA_2")
expColon<-read.table("./data_in/symbol.txt",header = T,check.names = F,sep="\t",stringsAsFactors = F)
clinColon<-read.table("./data_in/clinical.txt",header = T,check.names = F,sep="\t",stringsAsFactors = F)
expColon<-expColon[!grepl("\\.",expColon$ID),]
# reomove genes with counts less than 10 in more than 90% samples
expColon<-expColon[rowSums(expColon<10)<ncol(expColon)*0.9, ]
dim(expColon)
expColon=as.matrix(expColon)
rownames(expColon)=expColon[,1]
exp=expColon[,2:ncol(expColon)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
library(limma)
data=avereps(data)
data<-as.data.frame(data)
colnames(data)<-gsub("\\.","-",colnames(data))
data_tumor<-data[,grepl("^0",sapply(strsplit(colnames(data),"-"),'[',4))]
colnames(data_tumor)<-gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",colnames(data_tumor))
data_tumorUniq<-data_tumor[,!duplicated(colnames(data_tumor), fromLast=T)]
sampleOverlap<-intersect(colnames(data_tumorUniq),as.character(clinColon$ID))
data_tumorUniq_clinCommon<-data_tumorUniq[,sampleOverlap]
row.names(clinColon)<-clinColon$ID
clinColon_expCommon<-clinColon[sampleOverlap,]
write.table(data_tumorUniq_clinCommon,"./data_out/symbol_uniqueGene_tumorUniq_clinCommon.txt",row.names = T,quote = F,sep="\t")
write.table(clinColon_expCommon,"./data_out/clinColon_expCommon.txt",row.names = F,quote = F,sep="\t")
## ================
## WGCNA
## ================
### loading parametre
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()
exprMat <- data_tumorUniq_clinCommon
type = "unsigned"
corType = "non-pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
## Loading data
exprMat <- "./data_out/symbol_uniqueGene_tumorUniq_clinCommon.txt"
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T,
quote="", comment="", check.names=F)
dim(dataExpr)
dataExpr[c(1:6),c(1:6)]
## Data selection
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad >
max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr <- as.data.frame(t(dataExprVar))
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:",
paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:",
paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
# Remove the offending genes and samples from the data:
dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
## soft cutoff selection
sampleTree = hclust(dist(dataExpr), method = "average")
# Fig
# wgcna1_samples_hclust.pdf
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",cex=.2)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
networkType=type, verbose=5)
clust = cutreeStatic(sampleTree, cutHeight = 50000, minSize = 10)
table(clust)
row.names(dataExpr[clust==0, ])
keepSamples = (clust==1)
dataExpr = dataExpr[keepSamples, ]
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
# Re-cluster samples
sampleTree2 = hclust(dist(dataExpr), method = "average")
# wgcna1.2_samples_hclust.pdf
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="",cex=.2)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
networkType=type, verbose=5)
# wgcna2_power_sft.pdf
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",
ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
cex=cex1, col="red")
power = sft$powerEstimate
power
## network construction
nGenes
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
TOMType = type, minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs=TRUE, corType = corType,
maxPOutliers=maxPOutliers, loadTOMs=TRUE,
saveTOMFileBase = paste0(exprMat, ".tom"),
verbose = 3)
table(net$colors)
## hierarchical clustering models
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# wgcna3_hierarchicalClusteringModels.pdf
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
"Module colors",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
## Correlation among modules
MEs = net$MEs
MEs_col = MEs
MEs
colnames(MEs_col) = paste0("ME", labels2colors(
as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# wgcna4_correlationModuleEigengene.pdf
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
marDendro = c(3,3,2,4),
marHeatmap = c(3,4,2,2), plotDendrograms = T,
xLabelsAngle = 90)
## export network for Cytoscape
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
cyt = exportNetworkToCytoscape(TOM,
edgeFile = paste(exprMat, ".edges.txt", sep=""),
nodeFile = paste(exprMat, ".nodes.txt", sep=""),
weighted = TRUE, threshold = 0,
nodeNames = probes, nodeAttr = moduleColors)
## link phenotype data
trait <- "data_out/clinColon_expCommon.txt"
if(trait != "") {
traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
check.names=FALSE, comment='',quote="")
sampleName = rownames(dataExpr)
traitData = traitData[match(sampleName, rownames(traitData)), ]
}
if (corType=="pearsoon") {
modTraitCor = cor(MEs_col, traitData, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
modTraitCor = modTraitCorP$bicor
modTraitP   = modTraitCorP$p
}
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

# wgcna5_phenotypeLinked.pdf
par(mar=c(5,6,4,1)+.1)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData),
yLabels = colnames(MEs_col),
cex.lab = 0.5,
ySymbols = colnames(MEs_col), colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix, setStdMargins = FALSE,
cex.text = 0.5, zlim = c(-1,1),
main = paste("Module-trait relationships"))
corType
## Correlate one phenotype(fustat) with one module(MElightyellow)
if (corType=="pearson") {
geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(
as.matrix(geneModuleMembership), nSamples))
} else {
geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
geneModuleMembership = geneModuleMembershipA$bicor
MMPvalue   = geneModuleMembershipA$p
}
if (corType=="pearson") {
geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
geneTraitP = as.data.frame(corPvalueStudent(
as.matrix(geneTraitCor), nSamples))
} else {
geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
geneTraitCor = as.data.frame(geneTraitCorA$bicor)
geneTraitP   = as.data.frame(geneTraitCorA$p)
}

ME_colors<-sub("ME","",names(MEs_col))
for (module in ME_colors) {
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
moduleGenes = moduleColors == module
write.table(as.data.frame(geneModuleMembership[moduleGenes, module_column]),paste("./data_out/ModulePhenoCorGenes_",module,".txt",sep=""),row.names = T,quote = F,sep="\t")
}
for (pheno in names(traitData)){
pheno_column = match(pheno,colnames(traitData))
pdf(paste("./figure/Correlation_",pheno,".pdf",sep=""))
for (module in ME_colors) {
modNames = substring(colnames(MEs_col), 3)
module_column = match(module, modNames)
moduleGenes = moduleColors == module
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
abs(geneTraitCor[moduleGenes, pheno_column]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("Gene significance for", pheno),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}
dev.off()
}
rowSums(abs(modTraitCor[,-4]))
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/pip4_step1_WGCNA.Rhistory")
