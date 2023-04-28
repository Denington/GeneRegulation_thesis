library(rtracklayer)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
setwd('~/Documents/Projects/Accelerated_Regions/')
#chrom.sizes <- read.table('/Users/jk2014/Downloads/Moving_stuff_to_laptop/Presentations_documents/hg19.chrom.sizes')
for(i in 3:12){
  x <- import.wig(paste0(i,'.accels.including.accs.fullcons.wig'))
  x <- x[!is.na(x$score)]
  export.wig(unique(x),paste0(i,'.accels.including.accs.fullcons.wig'),genome='hg19')
  wigToBigWig(paste0(i,'.accels.including.accs.fullcons.wig'),seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19),clip = T )

  x <- import.wig(paste0(i,'.accels.all_but.fullcons.wig'))
  x <- x[!is.na(x$score)]
  export.wig(unique(x),paste0(i,'.accels.all_but.fullcons.wig'),genome='hg19')
  wigToBigWig(paste0(i,'.accels.all_but.fullcons.wig'),seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19),clip = T )
  
  
}

for(i in 3:12){
  x <- import.wig(paste0(i,'.accels.including.accs.mammal.wig'))
  x <- x[!is.na(x$score)]
  export.wig(unique(x),paste0(i,'.accels.including.accs.mammal.wig'),genome='hg19')
  wigToBigWig(paste0(i,'.accels.including.accs.mammal.wig'),seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19),clip = T )
  
  x <- import.wig(paste0(i,'.accels.all_but.mammal.wig'))
  x <- x[!is.na(x$score)]
  export.wig(unique(x),paste0(i,'.accels.all_but.mammal.wig'),genome='hg19')
  wigToBigWig(paste0(i,'.accels.all_but.mammal.wig'),seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19),clip = T )
  
  
}

###
for(i in 2:1){
  xx<-fread(paste0(i,'.accels.including.accs.fullcons.bed'))
  xx.gr <- makeGRangesFromDataFrame(xx,seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)
  xx.gr$V4 <- NULL
  xx.gr <- xx.gr[!is.na(xx.gr$V5)]
  xx.gr$score <- xx.gr$V5
  xx.gr$V5 <- NULL
  end(xx.gr) <- end(xx.gr) -1
  export.wig(unique(xx.gr),paste0(i,'.accels.including.accs.fullcons.wig'),genome='hg19')
  wigToBigWig(paste0(i,'.accels.including.accs.fullcons.wig'),seqinfo = seqinfo(BSgenome.Hsapiens.UCSC.hg19),clip = T)
}