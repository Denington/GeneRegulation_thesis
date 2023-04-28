#
setwd('~/Documents/Projects/Accelerated_Regions/')
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)


hg.m <- makeGRangesFromDataFrame(read.table('Acceleration_DNase_work/human_mouse_pairedscores.bed',header = T),keep.extra.columns = T)
hg.m <- hg.m[unique(queryHits(findOverlaps(hg.m,cnes.gr)))]
hg.acc <- import.bed('generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/hg19.bed')
mm.acc <- import.bed('generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/mm10.bed')
hg.bgc <- import.bed('generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/hg19.GCbias.bed')
mm.bgc <- import.bed('generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/mm10.GCbias.bed')

hg.m$accel <- 'None'

hg.m[unique(queryHits(findOverlaps(hg.m,mm.acc)))]$accel <- 'Mouse'
hg.m[unique(queryHits(findOverlaps(hg.m,hg.acc)))]$accel <- 'Human'

ggplot(as.data.frame(hg.m),aes(human_cell_lines,color=accel)) + stat_ecdf(geom  = 'line') +
  theme_classic() + ylab('Proportion of elements') + xlab('Number of human cell lines elements are active in') +
  ggsci::scale_color_jco() + theme(text = element_text(size = 20))

ggplot(as.data.frame(hg.m),aes(mouse_cell_lines,color=accel)) + stat_ecdf(geom  = 'line') +
  theme_classic() + ylab('Proportion of elements') + xlab('Number of mouse cell lines elements are active in') +
  ggsci::scale_color_jco() + theme(text = element_text(size = 20))

# accelerated regions are significantly more tissue specific in the species tested, but not for other species
# and note that human accels are no more likely to be 'off'-no longer active anywhere
#median and lower quartile numbers fairly similar- differences are less likely to be active everywhere
#whilst for mouse accels, see a higher proportion 'switched off'
#also ars ars aren't as tissue specific in other species

hg.m$gBGC <- 'None'

hg.m[unique(queryHits(findOverlaps(hg.m,mm.bgc)))]$gBGC <- 'Mouse'
hg.m[unique(queryHits(findOverlaps(hg.m,hg.bgc)))]$gBGC <- 'Human'

ggplot(as.data.frame(hg.m),aes(human_cell_lines,color=gBGC)) + stat_ecdf(geom  = 'line') +
  theme_classic() + ylab('Proportion of elements') + xlab('Number of human cell lines elements are active in') +
  ggsci::scale_color_jco()  + theme(text = element_text(size = 20))
ggplot(as.data.frame(hg.m),aes(mouse_cell_lines,color=gBGC)) + stat_ecdf(geom  = 'line') +
  theme_classic() + ylab('Proportion of elements') + xlab('Number of mouse cell lines elements are active in') +
  ggsci::scale_color_jco() + theme(text = element_text(size = 20))

# gc bias regions, if anything, are more active

cnes <- import.bed('CNEs.July17.bed')
accel_beds <- Sys.glob('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/*bed')
#accel_beds <- Sys.glob('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/*bed')

# need more species if using the nongc bias here
acc <- as.data.frame(cnes)[,c(1,2,3)]
for(i in accel_beds){
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  acc[,sp] <- 0
  accel.bed <- import.bed(i)
  #
  
  acc[unique(queryHits(findOverlaps(cnes,accel.bed))),][,sp] <- 1
  
}
rs <- rowSums(acc[,-c(1,2,3)])
acc.counts <- cnes
acc.counts$n <- rs

hg.m$RAIML_number <- 'None'

hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n == 1])))]$RAIML_number <- '1'
hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n == 2])))]$RAIML_number <- '2'
hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n == 3])))]$RAIML_number <- '3'
hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n == 4])))]$RAIML_number <- '4'
hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n > 4])))]$RAIML_number <- '5+'

ggplot(as.data.frame(hg.m),aes(human_cell_lines,color=RAIML_number)) + stat_ecdf(geom  = 'line') + theme_classic() +
  ylab('Proportion of elements') + xlab('Number of human cell lines elements are active in') + 
  ggsci::scale_color_jco() + theme(text = element_text(size = 20))
ggplot(as.data.frame(hg.m),aes(mouse_cell_lines,color=RAIML_number)) + stat_ecdf(geom  = 'line') + theme_classic() + 
  ylab('Proportion of elements') + xlab('Number of mouse cell lines elements are active in') +
  ggsci::scale_color_jco() + theme(text = element_text(size = 20))


#same,bgc
accel_beds <- Sys.glob('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/*bed')
#accel_beds <- Sys.glob('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/*bed')

acc <- as.data.frame(cnes)[,c(1,2,3)]
for(i in accel_beds){
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  acc[,sp] <- 0
  accel.bed <- import.bed(i)
  #
  
  acc[unique(queryHits(findOverlaps(cnes,accel.bed))),][,sp] <- 1
  
}
rs <- rowSums(acc[,-c(1,2,3)])
acc.counts <- cnes
acc.counts$n <- rs

hg.m$gBGC_number <- 'None'

hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n == 1])))]$gBGC_number <- '1'
hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n == 2])))]$gBGC_number <- '2'
hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n == 3])))]$gBGC_number <- '3'
hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n >= 4])))]$gBGC_number <- '4+'
#hg.m[unique(queryHits(findOverlaps(hg.m,acc.counts[acc.counts$n > 4])))]$gBGC_number <- '5+'

ggplot(as.data.frame(hg.m),aes(human_cell_lines,color=gBGC_number)) + stat_ecdf(geom  = 'line') + theme_classic() +
  ylab('Proportion of elements') + xlab('Number of human cell lines elements are active in') + 
  ggsci::scale_color_jco() + theme(text = element_text(size = 20))
ggplot(as.data.frame(hg.m),aes(mouse_cell_lines,color=gBGC_number)) + stat_ecdf(geom  = 'line') + theme_classic() + 
  ylab('Proportion of elements') + xlab('Number of mouse cell lines elements are active in') +
  ggsci::scale_color_jco() + theme(text = element_text(size = 20))

# regions with gBGC multiaccs are *less* tissue specific?? or at least not very different


#see if reallllly accelerated mouse regions are deleterious

h.Acc <-makeGRangesFromDataFrame(read.table('generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/hg19.bed'),seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)
mm.Acc <-makeGRangesFromDataFrame(read.table('generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/mm10.bed'),seqnames.field = 'V1',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)

mm.Acc$bins <- cut(mm.Acc$V4, breaks=c(0,4.4,5.35,6.8,9.38,15.816,100), labels=c("0%", "20%","40%","60%",'80%','95%+'))
h.Acc$bins <- cut(h.Acc$V4, breaks=c(0,6.89,7.46,8.14,9.5,12.606,100), labels=c("0%", "20%","40%","60%",'80%','95%+'))

quant.lrts <- function(a,cls){
  x0 <- cls[unique(queryHits(findOverlaps(cls,a[a$bins == '0%'])))] ; x0$quantile <- '0-20%'
  x2<- cls[unique(queryHits(findOverlaps(cls,a[a$bins == '20%'])))] ; x2$quantile <- '20-40%'
  x4 <- cls[unique(queryHits(findOverlaps(cls,a[a$bins == '40%'])))] ; x4$quantile <- '40-60%'
  x6 <- cls[unique(queryHits(findOverlaps(cls,a[a$bins == '60%'])))] ; x6$quantile <- '60-80%'
  x8 <- cls[unique(queryHits(findOverlaps(cls,a[a$bins == '80%'])))] ; x8$quantile <- '80-95%'
  x95 <- cls[unique(queryHits(findOverlaps(cls,a[a$bins == '95%+'])))] ; x95$quantile <- '95%+'
  return(c(x0,x2,x4,x6,x8,x95))
  
}

mouseLRTs <- as.data.frame(quant.lrts(mm.Acc,hg.m))

ggplot(mouseLRTs,aes(quantile,mouse_cell_lines)) + geom_boxplot() + theme_classic() + coord_cartesian(ylim=c(0,20))
ggplot(mouseLRTs,aes(quantile,human_cell_lines)) + geom_boxplot() + theme_classic() + coord_cartesian(ylim=c(0,20))
ggplot(mouseLRTs,aes(mouse_cell_lines,color=quantile)) + stat_ecdf(geom = 'line') + theme_bw() + 
  ylab('Proportion of elements') + xlab('Number of mouse cell lines elements are active in') + ggsci::scale_color_jco()
ggplot(mouseLRTs,aes(human_cell_lines,color=quantile)) + stat_ecdf(geom = 'line') + theme_bw() + 
  ylab('Proportion of elements') + xlab('Number of human cell lines elements are active in') +ggsci::scale_color_jco()

humanLRTs <- as.data.frame(quant.lrts(h.Acc,hg.m))

ggplot(humanLRTs,aes(mouse_cell_lines,color=quantile)) + stat_ecdf(geom = 'line') + theme_bw() + 
  ylab('Proportion of elements') + xlab('Number of mouse cell lines elements are active in') + ggsci::scale_color_jco()
ggplot(humanLRTs,aes(human_cell_lines,color=quantile)) + stat_ecdf(geom = 'line') + theme_bw() + 
  ylab('Proportion of elements') + xlab('Number of human cell lines elements are active in') +ggsci::scale_color_jco()


#
zero.actv <-  mm.Acc[unique(queryHits(findOverlaps(mm.Acc, hg.m[hg.m$mouse_cell_lines == 0])))]
low.actv <-  mm.Acc[unique(queryHits(findOverlaps(mm.Acc, hg.m[hg.m$mouse_cell_lines > 0 & hg.m$mouse_cell_lines < 5])))]
some.actv <- mm.Acc[unique(queryHits(findOverlaps(mm.Acc, hg.m[hg.m$mouse_cell_lines > 5 & hg.m$mouse_cell_lines < 15])))]
lots.actv <- mm.Acc[unique(queryHits(findOverlaps(mm.Acc, hg.m[hg.m$mouse_cell_lines > 15])))]

mm.lrt.v.c <- rbind( data.frame(Activity=rep('None',length(zero.actv$V4)),LRT=zero.actv$V4),
       data.frame(Activity=rep('Low',length(low.actv$V4)),LRT=low.actv$V4),
       data.frame(Activity=rep('Medium',length(some.actv$V4)),LRT=some.actv$V4),
       data.frame(Activity=rep('High',length(lots.actv$V4)),LRT=lots.actv$V4)
)

mm.lrt.v.c$Activity <- factor(mm.lrt.v.c$Activity,levels = c('None','Low','Medium','High'))
 
ggplot(mm.lrt.v.c,aes(LRT,color=Activity)) + stat_ecdf() + scale_x_continuous(trans = 'log10') +
  theme_classic() + ylab('Proportion of elements') + xlab('Number of mouse cell lines elements are active in') +
  ggsci::scale_color_jco()


ggplot(mm.lrt.v.c,aes(Activity,LRT,fill=Activity)) + geom_boxplot(outlier.alpha = 0.1) + scale_y_continuous(trans = 'log10') +
  theme_classic() + ggsci::scale_fill_jco() + ggpubr::stat_compare_means(comparisons = list(c('None','High'))) +
  theme(text = element_text(size = 15))
#flip box plots around so it's like proportion that are 0 cell lines active... etc
zero.actv <- hg.m[hg.m$mouse_cell_lines == 0]
some.actv <- hg.m[hg.m$mouse_cell_lines > 0]

mza <- length(zero.actv[unique(queryHits(findOverlaps(zero.actv,mm.acc)))])
msa <- length(some.actv[unique(queryHits(findOverlaps(some.actv,mm.acc)))])
cza <- length(zero.actv[unique(queryHits(findOverlaps(zero.actv,cnes.gr)))]) - mza
csa <- length(some.actv[unique(queryHits(findOverlaps(some.actv,cnes.gr)))]) - msa

fisher.test(matrix(c(mza,msa,cza,csa),ncol = 2))
#
zero.actv <- hg.m[hg.m$human_cell_lines == 0]
some.actv <- hg.m[hg.m$human_cell_lines > 0]

mza <- length(zero.actv[unique(queryHits(findOverlaps(zero.actv,mm.acc)))])
msa <- length(some.actv[unique(queryHits(findOverlaps(some.actv,mm.acc)))])
cza <- length(zero.actv[unique(queryHits(findOverlaps(zero.actv,cnes.gr)))]) - mza
csa <- length(some.actv[unique(queryHits(findOverlaps(some.actv,cnes.gr)))]) - msa

fisher.test(matrix(c(mza,msa,cza,csa),ncol = 2))
#
zero.actv <- hg.m[hg.m$human_cell_lines == 0]
some.actv <- hg.m[hg.m$human_cell_lines > 0]

mza <- length(zero.actv[unique(queryHits(findOverlaps(zero.actv,hg.acc)))])
msa <- length(some.actv[unique(queryHits(findOverlaps(some.actv,hg.acc)))])
cza <- length(zero.actv[unique(queryHits(findOverlaps(zero.actv,cnes.gr)))]) - mza
csa <- length(some.actv[unique(queryHits(findOverlaps(some.actv,cnes.gr)))]) - msa

fisher.test(matrix(c(mza,msa,cza,csa),ncol = 2))
#
#

######
#########
mycomps <- list(c('HAR','CNE'),c('CNE','MAR'))
allt <- data.frame(d= human.dnase$`spinal.cord-fetal`, sp='CNE')
hgtf <- data.frame(d= human.dnase[unique(queryHits(findOverlaps(human.dnase,hg.acc)))]$`spinal.cord-fetal`,sp='HAR')
mmtf <- data.frame(d= human.dnase[unique(queryHits(findOverlaps(human.dnase,mm.acc)))]$`spinal.cord-fetal`,sp='MAR')
ggplot(rbind(allt,hgtf,mmtf),aes(sp,d,fill=sp)) + geom_violin() + ylab('DNase signal (all elements)') + xlab('Type of element')+
  ggpubr::stat_compare_means(comparisons = mycomps) + theme_classic() + ggsci::scale_fill_jco()
#ggplot(rbind(allt,hgtf,mmtf),aes(d,color=sp)) + stat_ecdf() + theme_bw() + ggsci::scale_color_jco()

ggplot(rbind(allt[allt$d !=0,],hgtf[hgtf$d !=0,],mmtf[mmtf$d !=0,]),aes(sp,d,fill=sp)) + geom_violin() + 
  ggpubr::stat_compare_means(comparisons = mycomps) + theme_classic() + ggsci::scale_fill_jco() + 
  ylab('DNase signal (active elements)') + xlab('Type of element')

#ggplot(rbind(allt[allt$d !=0,],hgtf[hgtf$d !=0,],mmtf[mmtf$d !=0,]),aes(d,color=sp)) + stat_ecdf() + theme_bw() + ggsci::scale_color_jco()

allt <- data.frame(d= mouse.dnase$`retina-embryonic`, sp='CNE')
hgtf <- data.frame(d= mouse.dnase[unique(queryHits(findOverlaps(mouse.dnase,hg.acc)))]$`retina-embryonic`,sp='HAR')
mmtf <- data.frame(d= mouse.dnase[unique(queryHits(findOverlaps(mouse.dnase,mm.acc)))]$`retina-embryonic`,sp='MAR')
ggplot(rbind(allt,hgtf,mmtf),aes(sp,d,fill=sp)) + geom_violin() + ylab('DNase signal (all elements)') + xlab('Type of element')+
  ggpubr::stat_compare_means(comparisons = mycomps) + theme_classic() + ggsci::scale_fill_jco()
#ggplot(rbind(allt,hgtf,mmtf),aes(d,color=sp)) + stat_ecdf() + theme_bw() + ggsci::scale_color_jco()

ggplot(rbind(allt[allt$d !=0,],hgtf[hgtf$d !=0,],mmtf[mmtf$d !=0,]),aes(sp,d,fill=sp)) + geom_violin() + 
  ggpubr::stat_compare_means(comparisons = mycomps) + theme_classic() + ggsci::scale_fill_jco() + 
  ylab('DNase signal (active elements)') + xlab('Type of element')
#ggplot(rbind(allt[allt$d !=0,],hgtf[hgtf$d !=0,],mmtf[mmtf$d !=0,]),aes(d,color=sp)) + stat_ecdf() + theme_bw() + ggsci::scale_color_jco()



######## Addendum: sampling tissues from human data down to 64 cell types (like mouse)
######## and seeing if get the same increase in 0-cell types active for HARs

hg_dnase.m = read.table('~/Documents/Projects/Accelerated_Regions/Acceleration_DNase_work/human+0.5rep.collated.matrix.any.bed',sep='\t',header = T)
colLen <- length(colnames(hg_dnase.m))

hgTest <- hg_dnase.m[,4:colLen]
hgTest[hgTest > 1] <- 1
hgTest[hgTest < 1] <- 0

colNam <- colnames(hgTest)
keepCols <- sample(colnames(hgTest),64,replace = F)
hgTest <- hgTest[,keepCols]

sampled_cell_counts <- rowSums(hgTest)

hg_dnase_ss <- hg_dnase.m[,c('chrom','start','end')]
hg_dnase_ss$human_cell_lines <- sampled_cell_counts

hg_dnase_ss.gr <- makeGRangesFromDataFrame(hg_dnase_ss,seqnames.field = 'chrom',keep.extra.columns = T)
hg_dnase_ss.gr$accel <- 'None'

hg_dnase_ss.gr[unique(queryHits(findOverlaps(hg_dnase_ss.gr,mm.acc)))]$accel <- 'Mouse'
hg_dnase_ss.gr[unique(queryHits(findOverlaps(hg_dnase_ss.gr,hg.acc)))]$accel <- 'Human'

ggplot(as.data.frame(hg_dnase_ss.gr),aes(human_cell_lines,color=accel)) + stat_ecdf(geom  = 'line') +
  theme_classic() + ylab('Proportion of elements') + xlab('Number of human cell lines elements are active in') +
  ggsci::scale_color_jco() + theme(text = element_text(size = 10)) 


####
zero.actv <- hg_dnase_ss.gr[hg_dnase_ss.gr$human_cell_lines == 0]
some.actv <- hg_dnase_ss.gr[hg_dnase_ss.gr$human_cell_lines > 0]

mza <- length(zero.actv[unique(queryHits(findOverlaps(zero.actv,mm.acc)))])
msa <- length(some.actv[unique(queryHits(findOverlaps(some.actv,mm.acc)))])
cza <- length(zero.actv[unique(queryHits(findOverlaps(zero.actv,cnes.gr)))]) - mza
csa <- length(some.actv[unique(queryHits(findOverlaps(some.actv,cnes.gr)))]) - msa

fisher.test(matrix(c(mza,msa,cza,csa),ncol = 2))


