library("ATACseqQC")
library("Rsamtools")
library("ChIPpeakAnno")
library("TxDb.Dmelanogaster.UCSC.dm6.ensGene")
library("BiocIO")

bamfiles = c("A4_ED_rep2.sorted.chrX.bam", "A7_WD_rep1.sorted.chrX.bam", 
             "A4_ED_rep3.sorted.chrX.bam", "A7_WD_rep2.sorted.chrX.bam", 
             "A4_ED_rep4.sorted.chrX.bam", "A7_WD_rep3.sorted.chrX.bam")

outPath = "splited"
dir.create(outPath, showWarnings = FALSE)

for (bamfile in bamfiles) {

  bamfile.labels = gsub(".bam", "", basename(bamfile))
  
  fragSize = fragSizeDist(bamfile, bamfile.labels)
  
  estimateLibComplexity(readsDupFreq(bamfile, index=bamfile))
  
  tiff(paste0("figure_", bamfile.labels, ".tiff"), width = 7, height = 6, units = "in", res = 600)
  fragSize = fragSizeDist(bamfile, bamfile.labels)
  graphics.off()

  gal = readBamFile(bamfile, asMates=TRUE, bigFile=TRUE)
  gal1 = shiftGAlignmentsList(gal)
  
  shiftedBamfile = file.path(outPath, paste0(bamfile.labels, "_shifted.bam"))
  export(gal1, shiftedBamfile)

  txdbmm = makeTxDbFromGFF(file = "/pub/nshiroon/EE283/DATA/ATACseq/ref/dm6.ensGene.gtf", format = 'gtf')
  txs = transcripts(txdbmm)
  TSS = promoters(txs, upstream=0, downstream=1)
  TSS = unique(TSS)
  
  objs = splitGAlignmentsByCut(gal1, txs = txs, outPath = outPath)
  
  null = writeListOfGAlignments(objs, outPath)
  
  bamfiles_split = file.path(outPath, c("NucleosomeFree.bam", "mononucleosome.bam", 
                                        "dinucleosome.bam", "trinucleosome.bam"))
  
  librarySize = estLibSize(bamfiles_split)
  
  sigs = enrichedFragments(bamfiles_split, TSS = TSS, librarySize = librarySize, 
                           seqlev = "chrX", TSS.filter = 0.5, n.tile = 101, 
                           upstream = 1010, downstream = 1010)
  
  sigs.log2 <- lapply(sigs, function(.ele) log2(.ele + 1))
  
  featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width = 1010 + 1010),
                        zeroAt = 0.5, n.tile = 101)
}

