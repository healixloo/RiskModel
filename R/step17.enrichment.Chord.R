
#!/usr/bin/env Rscript
setwd("/Users/jlu/Desktop/Pro_Mine/Paper_TCGA_2/data_out")
library(circlize)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

david_enrichment<-read.csv(args[1],head=T,sep="\t",check.names = F)
david_enrichment<-david_enrichment[david_enrichment$FDR<0.05,]

david_enrichment_long<-tidyr::separate_rows(david_enrichment,Genes,sep=",")
david_enrichment_mat<-table(david_enrichment_long[,c("Term","Genes")])
rownames(david_enrichment_mat)<-str_split_fixed(rownames(david_enrichment_mat),"~",2)[,2]

pdf(paste(args[1],"pdf",sep="."),width=10,height=10)
grid.col <- setNames(rainbow(length(unlist(dimnames(david_enrichment_mat)))), union(rownames(david_enrichment_mat), colnames(david_enrichment_mat)))
par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))
# now, the image with rotated labels
chordDiagram(david_enrichment_mat, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
xlim = get.cell.meta.data("xlim")
ylim = get.cell.meta.data("ylim")
sector.name = get.cell.meta.data("sector.index")
circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
circos.axis(h = "top", labels.cex = 0.02, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
dev.off()

