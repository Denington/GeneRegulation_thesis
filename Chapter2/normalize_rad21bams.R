# fiddling with ways to normalize the rad21 bams, which have 2-fold issuea
# 1- comparing WT1 vs CTCF-/- normalize on peaks? But expect a bunch of peaks to be lost
# 2- WT4 (both replicates) have significantly higher signal to noise ratio
# so probably better to normalize at peaks
#could try having a WT vs CTCF, and a WT1 vs WT4 w different methods?
# will note at bottom of script downstream commands


wt1.peaks <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/Rad21_ChIPseq_in_CD69negDP/WT_CD69negDP/Rad21CD69negDPWTR1_mapped_sorted_RemoveDuplicates_peaks.subpeaks.bed',sep='\t',header = F),
                                          seqnames.field = 'V1',start.field='V2',end.field='V3')  
ctcfko.peaks <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/Rad21_ChIPseq_in_CD69negDP/CTCFKO_CD69negDP/Rad21CD69negDPCTCFKOR2_peaks.narrowPeak.bed',sep='\t',header = F),
                                         seqnames.field = 'V1',start.field='V2',end.field='V3')              

wt4.peaks <- makeGRangesFromDataFrame(read.table('~/Documents/Projects/Thymocyte_HiC/Thesis/Ya_chips/6_Rad21CD69pos4SPWTR2_mapped_sorted_RemoveDuplicates_MACS_bedGraph/Rad21CD69pos4SPWTR2_mapped_sorted_RemoveDuplicates_peaks.subpeaks.bed',sep='\t',header = T),
                                      seqnames.field = 'Chromosome',start.field='Start',end.field='End')

test.regions <- Reduce(subsetByOverlaps, list(wt1.peaks, wt4.peaks, ctcfko.peaks))
#Reduce(intersect, list(gr, gr1, gr2))

reg.counts <- assay(regionCounts(bam.files, test.regions, ext=frag.len, param=param)) 
#raw.counts <-  assay(filtered.data) #filtered.data@assays@data$counts

norm.factor   <- calcNormFactors(object = reg.counts, method = c("TMM"))
lib.size      <- colSums(reg.counts)
final.factor  <- norm.factor * lib.size

## as one typically scales reads "per million" we divide this factor by 10^6
## we also need the reciprocal value when using bamCoverage later, you'll see why,
## see *comment below
perMillion.factor <- (final.factor / 1000000)^-1

write.table(x = data.frame(Sample = bam.files,   NormFactor = perMillion.factor),
            file = "/Users/jk2014/Documents/Projects/Thymocyte_HiC/chip_seqs/bams/normFactor_normedBams/dif.rad21.normFactors.txt", sep="\t", quote = FALSE,
            col.names = TRUE, row.names = FALSE)

'''
for i in *.bam
  do
  
  ## extract norm.factor from the file we created in R for the current sample:
  NormFactor=$(grep "${i}" dif.rad21.normFactors.txt | cut -f2)
  
  ## use bamCoverage:
  bamCoverage -b $i -o ${i%.bam}_TMM.bigwig -bs 1 --scaleFactor ${NormFactor}
  done
'''


