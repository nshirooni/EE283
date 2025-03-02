library(tidyverse)
mytab = read_tsv("RNAseq384_SampleCoding.txt")
mytab
mytab2 <- mytab %>%
select(RILcode, TissueCode, Replicate, FullSampleName)
table(mytab2$RILcode)
table(mytab2$TissueCode)
table(mytab2$Replicate)

mytab2 <- mytab %>%
select(RILcode, TissueCode, Replicate, FullSampleName) %>%
filter(RILcode %in% c(21029, 21148, 21286, 21293, 21297, 22031,
22052, 22162, 22378, 22390)) %>%
filter(TissueCode %in% c("B", "E")) %>%
filter(Replicate == "0")

mytab2 <- mytab2 %>%
  mutate(FullSampleName = paste0(RILcode, "_", TissueCode, "_", Replicate))

for(i in 1:nrow(mytab2)){
cat("RNAseq/bam/",mytab2$FullSampleName[i],".bam\n",file="shortRNAseq.names
.txt",append=TRUE,sep='')
}
write_tsv(mytab2,"shortRNAseq.txt")
