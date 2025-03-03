library(DESeq2)
library(gplots)
library(RColorBrewer)
library(genefilter)

sampleInfo = read.table("shortRNAseq.txt", header = TRUE)
sampleInfo$FullSampleName = as.character(sampleInfo$FullSampleName)
countdata = read.table("fly_counts.txt", header=TRUE, row.names=1)
countdata = countdata[ ,6:ncol(countdata)]
temp = colnames(countdata)
temp = gsub("RNAseq.bam.","",temp)
temp = gsub(".bam","",temp)
colnames(countdata) = temp
cbind(temp,sampleInfo$FullSampleName,temp == sampleInfo$FullSampleName)
dds = DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo,
                             design=~TissueCode)
dds <- DESeq(dds)
res <- results( dds )
res
plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds )
hist( res$pvalue, breaks=20, col="grey" )
rld = rlog( dds )
each sample ID, and then roll your own analysis
head( assay(rld) )
mydata = assay(rld)
sampleDists = dist( t( assay(rld) ) )
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = rld$TissueCode
colnames(sampleDistMatrix) = NULL
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
print( plotPCA( rld, intgroup = "TissueCode") )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none",
           dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

library(ggplot2)
pval = res$pvalue
volcano_data = data.frame(
  log2fc = res$log2FoldChange,
  pval = res$pvalue,
  label = ifelse(pval < 0.05 & abs(res$log2FoldChange) > 1, rownames(res), "")
)

volcano_data$logpval = -log10(volcano_data$pval)

ggplot(volcano_data, aes(x=log2fc, y=logpval)) +
  geom_point(aes(color = (pval < 0.05 & abs(log2fc) > 1)), size=1) +
  scale_color_manual(values=c("black", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="blue") + 
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="blue") + 
  theme_minimal() +
  xlim(-4, 4) + 
  ylim(0, 20) + 
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 P-value") +
  theme(legend.position="none")
