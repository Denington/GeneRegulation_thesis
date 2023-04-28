##coverage heatmap plots of accels

library(genomation)
library(rtracklayer)
library(RColorBrewer)
library(ggsci)
library(data.table)
library(tidyverse)
library(magrittr)
setwd('~/Documents/Projects/Accelerated_Regions/')

#x <- read.table('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/grb_boundaries_genes.hg19.txt',col.names = c('seqnames','start','end','target_gene'))
x <- read.table('~/Documents/Projects/Accelerated_Regions/hg19.canFam3.98.50.processed.threshold.grbs.bed',skip = 1) # col.names = c('seqnames','start','end','target_gene'))

grbs.gr <- unique(makeGRangesFromDataFrame(x,keep.extra.columns = T,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3'))
grbs.gr <- grbs.gr[width(grbs.gr) > 250000]

#grbs.gr <- grbs.gr[order(as.factor(grbs.gr$target_gene))]
####
cnes <- import.bed('CNEs.July17.bed')
cnes <- unique(cnes)
human_accels <- import.bed('generalstuff/July Stuff/Accelerated_Regions.July2017/hg19.AcceleratedRegions.bed')
#
accels_per_cne <- fread("/Users/jk2014/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/perms grb stuff/grb enrichment/accelerated_table.July2017.tsv")
accels_per_cne <- accels_per_cne[,-c(4,5,10,11,15,16,19,20,21,24,25,26)]
accels_per_cne.g <- accels_per_cne %>% separate(coords,c("seqnames", "start"),c(':'))
accels_per_cne.g <- accels_per_cne.g %>% separate(start,c("start","end"),c('-'))
accels_per_cne.gr <- makeGRangesFromDataFrame(accels_per_cne.g,keep.extra.columns = T)

accels_per_cne <- accels_per_cne[unique(subjectHits(findOverlaps(accels_per_cne.gr,cnes))),]
#just filtering out some lower quality cnes
real_values <- accels_per_cne #[,c('coords','mm10','canFam3')]
#real_values <- real_values[,-c(3,4,6,7,11)]
real_counts<-rowSums(real_values[, 2:ncol(real_values)])

ars <- accels_per_cne.gr[real_counts > 0]
#
#try looking at all possible accelerated regions here
#all_accels.gr[rowSums(as.data.frame(all_accels.gr)[,-c(1,2,3,4,5)]) > 0,]
####

#too many bins I think
grb.Coverage.sm <- ScoreMatrixBin(coverage(grbs.gr),resize(grbs.gr[order(width(grbs.gr),decreasing = T)],width = 3000000,fix = 'center'),
                                  bin.num=500,strand.aware=F,type = 'max')
grb.Coverage.sm@.Data[which(grb.Coverage.sm@.Data>1)] = 1
heatMatrix(grb.Coverage.sm, xcoords = c(-1500000, 1500000),
           col = c('white','black'),main = 'GRB Coverage')
#### why is there a singly grb there's a double of?

CNE.Coverage.sm <- ScoreMatrixBin(coverage(cnes),resize(grbs.gr[order(width(grbs.gr),decreasing = T)],width = 3000000,fix = 'center'),
                                 bin.num=500,strand.aware=F,type = 'max')
#CNE.Coverage.sm@.Data <- round(CNE.Coverage.sm@.Data,0)
heatMatrix(CNE.Coverage.sm, xcoords = c(-1500000, 1500000),
           col = blues9,main = 'CNE Coverage at GRBs',winsorize = c(0,98))
#CNE.Coverage.sm@.Data[which(CNE.Coverage.sm@.Data> 0)] = 1
#M[which(M>0)] = 1
har.Coverage.sm <- ScoreMatrixBin(coverage(ars),resize(grbs.gr[order(width(grbs.gr),decreasing = T)],width = 3000000,fix = 'center'),
                                  bin.num=500,strand.aware=F,type='max')
#har.Coverage.sm@.Data <- round(har.Coverage.sm@.Data,0)
#har.Coverage.sm@.Data[which(har.Coverage.sm@.Data> 0)] = 1
#har.Coverage.sm@.Data <- har.Coverage.sm@.Data + CNE.Coverage.sm@.Data
#har.Coverage.sm@.Data[which(har.Coverage.sm@.Data== 3)] = 2
heatMatrix(har.Coverage.sm, xcoords = c(-1500000, 1500000),
           col = blues9,main = 'AR Coverage at GRBs',winsorize = c(0,98))


##### now for GC bias plots
