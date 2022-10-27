# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-06-ExportNetwork.pdf
setwd("~/Desktop/Pro_Mine/Paper_TCGA_2")
library(limma)
library(WGCNA)
library(reshape2)
library(stringr)
load("./data_in/wgcna_2.RData")
# Exporting to Cytoscape
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(dataExpr, power = 7)
# Select modules
modules = c("tan", "turquoise")
inModule = is.finite(match(moduleColors, modules))
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
genes = names(dataExpr)
modGenes = genes[inModule]
dimnames(modTOM) = list(modGenes, modGenes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""), nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""), weighted = TRUE,
threshold = 0.01, # 0.02
nodeNames = modGenes,
altNodeNames = modGenes,
nodeAttr = moduleColors[inModule])
savehistory("~/Desktop/Pro_Mine/Paper_TCGA_2/script/step1_WGCNA_add.R")
