library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(tidyverse)
setwd('~/Documents/Projects/Accelerated_Regions/')

cnes <- read.table('CNEs.July17.bed',sep='\t')
colnames(cnes) <- c('chrom','start','end')
cnes.gr <- import.bed('CNEs.July17.bed')

accel_beds <- Sys.glob('generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/*bed')
for(i in accel_beds){
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  accel.bed <- import.bed(i)
  
  cnes[,sp] <- 0
  cnes[,sp][unique(queryHits(findOverlaps(cnes.gr,accel.bed)))] <- 1
  
}

write.table(cnes,'no.bgc.accels.per.cne.tsv',sep='\t',quote = F,row.names = F)
#cnes <- cnes[,-c(9,11,15,16,17,18,20,21,22,25,27,28)]
#write.table(cnes,'14species.2022.accels.per.cne.tsv',sep='\t',quote = F,row.names = F)
