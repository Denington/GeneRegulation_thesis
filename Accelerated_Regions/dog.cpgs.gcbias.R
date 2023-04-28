library(rtracklayer)
library(Rsamtools)
library(genomation)
library(ggplot2)

setwd('~/Documents/Projects/Accelerated_Regions/DogBams.Sundman2020//')

#blacklist <- import.bed('../../Thesis/external_data/mm9_blacklist.bed')
#blacklist <- c(blacklist,bad.bins)

standard.chr <- paste0("chr", c(1:38))
#param <- readParam(minq=25, restrict=standard.chr)

bam.files <- c('134_Boxer_GBS_sorted_reheader.bam','134_Boxer_MeDIP_sorted_reheader.bam')
indexBam(bam.files)

###############
library(GenomicFeatures)
###

dogSession<-browserSession("UCSC")
genome(dogSession)<-"canFam3"

cpgIslands<-getTable(ucscTableQuery(dogSession, track="cpgIslandExt", table="cpgIslandExt"))
cpgIslands.gr<-GRanges(seqnames=cpgIslands$chrom, 
                       ranges=IRanges(start=cpgIslands$chromStart,
                                      end=cpgIslands$chromEnd),
                       gcNum=cpgIslands$gcNum,
                       cpgNum=cpgIslands$cpgNum)

##############

bgc <- import.bed('dog.gc.biasaccels.dogcoords.bed')
bgc <- reduce(bgc,ignore.strand = TRUE, min.gapwidth = 51L)
nobgc <- import.bed('dog.other.accels.dogcoords.bed')
nobgc <- reduce(nobgc,ignore.strand = TRUE, min.gapwidth = 51L)

bgc <- bgc[seqnames(bgc) %in% paste0('chr',1:38)]
nobgc <- nobgc[seqnames(nobgc) %in% paste0('chr',1:38)]
seqlevels(bgc) <- seqlevelsInUse(bgc)
seqlevels(nobgc) <- seqlevelsInUse(nobgc)
#
bgc <- bgc[width(bgc) > 90 & width(bgc) < 1000]
nobgc <- nobgc[width(nobgc) > 90 & width(nobgc) < 1000]

bgc.x <- resize(bgc,width(bgc)*20,fix = 'center')
nobgc.x <- resize(nobgc,width(nobgc)*20,fix = 'center')

bgc.sm <- ScoreMatrixBin(bam.files[[1]],bgc.x,bin.num = 30,strand.aware = F)
nobgc.sm <- ScoreMatrixBin(bam.files[[1]],nobgc.x,bin.num = 30,strand.aware = F)

plotMeta(ScoreMatrixList(List(bgc.sm,nobgc.sm)),winsorize = c(1,99), dispersion = 'se',  
         line.col = rainbow(20)[c(1,13)],dispersion.col = rainbow(20,alpha=0.35)[c(1,13)],
         xaxt="n",smoothfun=function(x) smooth.spline(x, spar=0.25),
         ylab='Average CpG methylation Coverage',xlab='',profile.names = c('Biased Gene Conversion','Other Accelerated Regions'))

axis(1, at=c(7.5,22.5), labels=c("Region Start", "Region End")) 
abline(v=c(7.5,22.5), lty=2)

bgc.df <- data.frame(vals = bgc.sm@.Data[,15],Type='Biased Gene Conversion')
no.df <- data.frame(vals = nobgc.sm@.Data[,15],Type='Other Dog accelerated Region')

ggplot(rbind(bgc.df,no.df),aes(Type,vals,fill=Type)) + geom_violin() + 
  coord_cartesian(ylim=c(0,15)) + ylab('CpG methylation') + theme_classic() +
  ggsci::scale_fill_jco() + ggpubr::stat_compare_means()

########
# CpG methylation levels are comprable at dog bgc and non bgc cpg islands
cpgI.gr <- cpgIslands.gr
x <- cpgI.gr[queryHits(findOverlaps(cpgI.gr,bgc))]
xx <- cpgI.gr[-queryHits(findOverlaps(cpgI.gr,bgc))]
xxx <- cpgI.gr[queryHits(findOverlaps(cpgI.gr,nobgc))]
seqlevels(x) <- seqlevelsInUse(x); xx <- keepSeqlevels(xx, paste0('chr',1:38), pruning.mode="coarse")
seqlevels(xx) <- seqlevelsInUse(xx); seqlevels(xxx) <- seqlevelsInUse(xxx)
cpgI.bgc.sm <- ScoreMatrixBin(bam.files[[1]],resize(x,4*width(x),fix = 'center'),bin.num = 30,strand.aware = F)
cpgI.no.sm <- ScoreMatrixBin(bam.files[[1]],resize(xx,4*width(xx),fix = 'center'),bin.num = 30,strand.aware = F)
cpgI.acc.sm <- ScoreMatrixBin(bam.files[[1]],resize(xxx,4*width(xxx),fix = 'center'),bin.num = 30,strand.aware = F)

plotMeta(ScoreMatrixList(list(cpgI.bgc.sm,cpgI.no.sm)),winsorize = c(2,98),
         profile.names = c('Biased Gene Conversion','Other CpG islands'),dispersion = 'se',
         line.col = rainbow(20)[c(1,13)],dispersion.col = rainbow(20,alpha=0.35)[c(1,13)],
         xaxt="n",smoothfun=function(x) smooth.spline(x, spar=0.25),
         ylab='Average CpG methylation Coverage',xlab='')

axis(1, at=c(7.5,22.5), labels=c("CpG Island Start", "CpG Island End")) 
abline(v=c(7.5,22.5), lty=2)


########################## same for mouse ################

bg <- import.bedGraph('~/Downloads/GSM1027571_DNA_CpG_coverage_E14_serum_LIF.bedGraph')
gbc <- import.bed('~/Downloads/mouse.gc.biasaccels.mm9coords.bed')
no.gbc <- import.bed('~/Downloads/mouse.other.accels.mm9coords.bed')
gbc <- gbc[order(width(gbc),decreasing = T)]
gbc <- gbc[width(gbc) < 1000 & width(gbc) > 90]

no.gbc <- no.gbc[order(width(no.gbc),decreasing = T)]
no.gbc <- no.gbc[width(no.gbc) < 1000 & width(no.gbc) > 90]

bgc.x <- resize(gbc,width(gbc)*10,fix = 'center')
nobgc.x <- resize(no.gbc,width(no.gbc)*10,fix = 'center')

bgc.sm <- ScoreMatrixBin(bg,bgc.x,bin.num = 30,strand.aware = F,weight.col = 'score')
nobgc.sm <- ScoreMatrixBin(bg,nobgc.x,bin.num = 30,strand.aware = F,weight.col = 'score')

plotMeta(ScoreMatrixList(List(bgc.sm,nobgc.sm)),winsorize = c(1,99), dispersion = 'se',  
         line.col = rainbow(20)[c(1,13)],dispersion.col = rainbow(20,alpha=0.35)[c(1,13)],
         xaxt="n",smoothfun=function(x) smooth.spline(x, spar=0.25),
         ylab='Average CpG methylation Coverage',xlab='',profile.names = c('Biased Gene Conversion','Other'))
axis(1, at=c(7.5,22.5), labels=c("Region Start", "Region End")) 
abline(v=c(7.5,22.5), lty=2)

bgc.df <- data.frame(vals = bgc.sm@.Data[,15],Type='Biased Gene Conversion')
no.df <- data.frame(vals = nobgc.sm@.Data[,15],Type='Other Mouse accelerated Region')

ggplot(rbind(bgc.df,no.df),aes(Type,vals,fill=Type)) + geom_violin() + 
  coord_cartesian(ylim=c(0,5)) + ylab('CpG methylation') + theme_classic() +
  ggsci::scale_fill_jco() + ggpubr::stat_compare_means()

#

genome(dogSession)<-"mm9"

cpgIslands<-getTable(ucscTableQuery(dogSession, track="cpgIslandExt", table="cpgIslandExt"))
cpgIslands.gr<-GRanges(seqnames=cpgIslands$chrom, 
                       ranges=IRanges(start=cpgIslands$chromStart,
                                      end=cpgIslands$chromEnd),
                       gcNum=cpgIslands$gcNum,
                       cpgNum=cpgIslands$cpgNum)

cpgI.gr <- cpgIslands.gr[seqnames(cpgIslands.gr) %in% paste0('chr',1:19)]
seqlevels(cpgI.gr) <- seqlevelsInUse(cpgI.gr)

cpgI.bgc.gr <- cpgI.gr[unique(queryHits(findOverlaps(cpgI.gr,gbc)))]
cpgI.nobgc.gr <- cpgI.gr[-unique(queryHits(findOverlaps(cpgI.gr,gbc)))]

cpgI.bgc.sm <- ScoreMatrixBin(bg,resize(cpgI.bgc.gr,4*width(cpgI.bgc.gr),fix = 'center'),bin.num = 30,strand.aware = F,weight.col = 'score')
cpgI.nobgc.sm <- ScoreMatrixBin(bg,resize(cpgI.nobgc.gr,4*width(cpgI.nobgc.gr),fix = 'center'),bin.num = 30,strand.aware = F,weight.col = 'score')
#plotMeta(cpgI.sm)


plotMeta(ScoreMatrixList(list(cpgI.bgc.sm,cpgI.nobgc.sm)),winsorize = c(2,98),
         profile.names = c('Biased Gene Conversion','No BGC'),dispersion = 'se',
         line.col = rainbow(20)[c(1,13)],dispersion.col = rainbow(20,alpha=0.35)[c(1,13)],
         xaxt="n",smoothfun=function(x) smooth.spline(x, spar=0.25),
         ylab='Average Coverage',xlab='')

axis(1, at=c(7.5,22.5), labels=c("CpG Island Start", "CpG Island End")) 
abline(v=c(7.5,22.5), lty=2)



# dog bgc regions are overwhelmingly promoter. And in other species? i think the dog thing is known tbf
#roller.dog <- read.table('~/Documents/Projects/Accelerated_Regions/Roller2021data/mergeRegRegions/Dog_regRegions_allTissue_parsed.txt')
#roller.dog$seqnames <- paste0('chr',roller.dog$V1)
#rol.dog.gr <- makeGRangesFromDataFrame(roller.dog,seqnames.field = 'seqnames',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)

#table(rol.dog.gr[unique(queryHits(findOverlaps(rol.dog.gr,nobgc)))]$V4)
#table(rol.dog.gr[unique(queryHits(findOverlaps(rol.dog.gr,bgc)))]$V4)

