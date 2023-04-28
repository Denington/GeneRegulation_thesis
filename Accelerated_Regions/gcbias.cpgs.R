#looking at how gcbias regions compare to cpg islands - 
#see that gcb accels are strongly enriched to be at cpg islands, and enriched to have already been cpg islands
# ie not just ...
# eg they're cpgislands in humans, which don't have gcb. Dogs etc known to have much higher density of CpGs, so maybe that's related
# gcb regions are also much closer to tss than non-gcb accels
# https://pubmed.ncbi.nlm.nih.gov/18477403/ fig1a  . https://pubmed.ncbi.nlm.nih.gov/19409480/ - dpg higher turnover at prom sequences also
# see all the above for mouse and human accels also, though nowhere near as strong  
#as only 7 mouse gcb regions overlap cpg islands
setwd('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/')
library(rtracklayer)
library(GenomicFeatures)
library(regioneR)
library(rtracklayer)
library(ggplot2)
###

cf.gcb <- import.bed('~/Downloads/dog.gc.biasaccels.dogcoords.bed')
cf.nogcb <- import.bed('~/Downloads/dog.other.accels.dogcoords.bed')

dogSession<-browserSession("UCSC")
genome(dogSession)<-"canFam3"

cpgIslands<-getTable(ucscTableQuery(dogSession, track="cpgIslandExt", table="cpgIslandExt"))
cpgIslands.gr<-GRanges(seqnames=cpgIslands$chrom, 
                         ranges=IRanges(start=cpgIslands$chromStart,
                                        end=cpgIslands$chromEnd),
                         gcNum=cpgIslands$gcNum,
                         cpgNum=cpgIslands$cpgNum)

doggenes <- getTable(ucscTableQuery(dogSession, track="NCBI RefSeq", table="ncbiRefSeq"))
doggenes.gr <- unique(GRanges(seqnames=doggenes$chrom, 
                       ranges=IRanges(start=doggenes$txStart,
                                      end=doggenes$txEnd),
                       name=doggenes$name2,
                       ))

plot(permTest(A= cf.gcb,B=cpgIslands.gr, universe=c(cf.gcb,cf.nogcb),ntimes = 500,
              randomize.function = resampleRegions,evaluate.function = numOverlaps, force.parallel = T))

m <- matrix(c(numOverlaps(cf.gcb,cpgIslands.gr,count.once = T), length(cf.gcb) - numOverlaps(cf.gcb,cpgIslands.gr,count.once = T),
numOverlaps(cf.nogcb,cpgIslands.gr,count.once = T), length(cf.nogcb) - numOverlaps(cf.nogcb,cpgIslands.gr,count.once = T)) ,ncol = 2)
fisher.test(m)

########dog CpG islands
summary(width(cpgIslands.gr))
summary(width(cpgIslands.gr[unique(queryHits(findOverlaps(cpgIslands.gr,cf.gcb)))]))
x1 <- width(cpgIslands.gr); x2 <- width(cpgIslands.gr[unique(queryHits(findOverlaps(cpgIslands.gr,cf.gcb)))])
zz<-rbind(data.frame(tt='non',wid=x1),data.frame(tt='bgc',wid=x2)) ;
ggplot(zz,aes(tt,wid,fill=tt)) + geom_boxplot() + scale_y_continuous(trans='log10')
########
#dog gcb accels are much more likely to be cpg islands

plot(permTest(A= cf.gcb,B=resize(doggenes.gr,10000,fix = 'center'), universe=c(cf.gcb,cf.nogcb),ntimes = 500,
              randomize.function = resampleRegions,evaluate.function = numOverlaps, force.parallel = T))

plot(permTest(A= cf.gcb,B=doggenes.gr, universe=c(cf.gcb,cf.nogcb),ntimes = 500,
              randomize.function = resampleRegions,evaluate.function = meanDistance, force.parallel = T))
#and are closer to genes

##################
cf.gcb.hgcoords <- import.bed('canFam3.GCbias.bed')
cf.nogcb.hgcoords <- import.bed('noGCbias/canFam3.bed')

humSession<-browserSession("UCSC")
genome(humSession)<-"hg19"

cpgIslands.hg <-getTable(ucscTableQuery(humSession, track="cpgIslandExt", table="cpgIslandExt"))
cpgIslands.hg.gr<-GRanges(seqnames=cpgIslands.hg$chrom, 
                       ranges=IRanges(start=cpgIslands.hg$chromStart,
                                      end=cpgIslands.hg$chromEnd),
                       gcNum=cpgIslands.hg$gcNum,
                       cpgNum=cpgIslands.hg$cpgNum)

hggenes <- getTable(ucscTableQuery(humSession, track="NCBI RefSeq", table="ncbiRefSeq"))
hggenes.gr <- unique(GRanges(seqnames=hggenes$chrom, 
                              ranges=IRanges(start=hggenes$txStart,
                                             end=hggenes$txEnd),
                              name=hggenes$name2,
))

plot(permTest(A= cf.gcb.hgcoords,B=cpgIslands.hg.gr, universe=c(cf.gcb.hgcoords,cf.nogcb.hgcoords),ntimes = 500,
              randomize.function = resampleRegions,evaluate.function = numOverlaps, force.parallel = T))

m <- matrix(c(numOverlaps(cf.gcb.hgcoords,cpgIslands.hg.gr,count.once = T), length(cf.gcb.hgcoords) - numOverlaps(cf.gcb.hgcoords,cpgIslands.hg.gr,count.once = T),
              numOverlaps(cf.nogcb.hgcoords,cpgIslands.hg.gr,count.once = T), length(cf.nogcb.hgcoords) - numOverlaps(cf.nogcb.hgcoords,cpgIslands.hg.gr,count.once = T)) ,ncol = 2)
fisher.test(m)

#dog gcb accels are significantly much more likely to be cpg islands even in human

plot(permTest(A= cf.gcb.hgcoords,B=resize(hggenes.gr,10000,fix = 'center'), universe=c(cf.gcb.hgcoords,cf.nogcb.hgcoords),ntimes = 500,
              randomize.function = resampleRegions,evaluate.function = numOverlaps, force.parallel = T))

plot(permTest(A= cf.gcb.hgcoords,B=hggenes.gr, universe=c(cf.gcb.hgcoords,cf.nogcb.hgcoords),ntimes = 500,
              randomize.function = resampleRegions,evaluate.function = meanDistance, force.parallel = T))



########
#getting CpG island density per megabase for a few species and comparing to number of (or proportion of ) BGC accelerations

get.CpG.density <- function(sp,bgc.accels,other.acc){
  Session<-browserSession("UCSC")
  genome(Session)<- sp
  
  cpgIslands <-getTable(ucscTableQuery(Session, track="cpgIslandExt", table="cpgIslandExt"))
  cpgIslands.gr<-GRanges(seqnames=cpgIslands$chrom, 
                            ranges=IRanges(start=cpgIslands$chromStart,
                                           end=cpgIslands$chromEnd),
                            gcNum=cpgIslands$gcNum,
                            cpgNum=cpgIslands$cpgNum)
  
  chrom.sizes <- read.table(paste0('~/Documents/Projects/Accelerated_Regions/chromSizes/',sp,'.chrom.sizes.txt'))
  norm.chroms <- length(chrom.sizes[chrom.sizes$V1 %in% paste0('chr',1:40),]$V1) + 2
  #
  cpgIslands.gr <- cpgIslands.gr[seqnames(cpgIslands.gr) %in% chrom.sizes$V1]
  chrom.sizes <- chrom.sizes[chrom.sizes$V1 %in% seqnames(cpgIslands.gr),]
  #
  gr.windows <- tileGenome(with(chrom.sizes, setNames(V2, V1)), tilewidth=100000,cut.last.tile.in.chrom=TRUE)
  gr.windows$val <- countOverlaps(gr.windows,cpgIslands.gr) 
  cpg.density <- mean(gr.windows$val)
  
  #totGenome <- sum(chrom.sizes$V2)
  #totCpG <- sum(width(cpgIslands.gr))
  #cova <- 100 * totCpG / totGenome
  bgc.acc <- length(import.bed(bgc.accels))
  other.acc <- length(import.bed(other.acc))
  #return(c(cova,bgc.acc/(bgc.acc+other.acc),sp))
  return(c(length(cpgIslands.gr),cpg.density,bgc.acc/(bgc.acc+other.acc),sp,norm.chroms))
}
#this bounces massively for dolphin depending on which method i use for coverage- because scaffolds aren't that big, tile of 1mb is too great
equ.cpg <- get.CpG.density('equCab3','equCab2.GCbias.bed','noGCbias/equCab2.bed')
hg.cpg <- get.CpG.density('hg19','hg19.GCbias.bed','noGCbias/hg19.bed')
panTro.cpg <- get.CpG.density('panTro4','panTro4.GCbias.bed','noGCbias/panTro4.bed')
rheMac3.cpg <- get.CpG.density('rheMac3','rheMac3.GCbias.bed','noGCbias/rheMac3.bed')
mm10.cpg <- get.CpG.density('mm10','mm10.GCbias.bed','noGCbias/mm10.bed')
rn5.cpg <- get.CpG.density('rn5','rn5.GCbias.bed','noGCbias/rn5.bed')
bosTau7.cpg <- get.CpG.density('bosTau7','bosTau7.GCbias.bed','noGCbias/bosTau7.bed')
turTru2.cpg <- get.CpG.density('turTru2','turTru2.GCbias.bed','noGCbias/turTru2.bed')
canFam3.cpg <- get.CpG.density('canFam3','canFam3.GCbias.bed','noGCbias/canFam3.bed')
felCat5.cpg <- get.CpG.density('felCat5','felCat5.GCbias.bed','noGCbias/felCat5.bed')
myoLuc2.cpg <- get.CpG.density('myoLuc2','myoLuc2.GCbias.bed','noGCbias/myoLuc2.bed')
susScr3.cpg <- get.CpG.density('susScr3','susScr3.GCbias.bed','noGCbias/susScr3.bed')
cavPor3.cpg <- get.CpG.density('cavPor3','cavPor3.GCbias.bed','noGCbias/cavPor3.bed')
#loxAfr3.cpg <- get.CpG.density('loxAfr3','loxAfr3.GCbias.bed','noGCbias/loxAfr3.bed')
#elephant has no cpg track

cpg.coverages <- as.data.frame(rbind(equ.cpg,hg.cpg,panTro.cpg,rheMac3.cpg,mm10.cpg,rn5.cpg,bosTau7.cpg,
                                     canFam3.cpg,felCat5.cpg,susScr3.cpg)) #cavPor3.cpg myoLuc2.cpg turTru2.cpg
cpg.coverages$V1 <- as.numeric(cpg.coverages$V1) ; cpg.coverages$V2 <- as.numeric(cpg.coverages$V2) ; 
cpg.coverages$V3 <- as.numeric(cpg.coverages$V3) ; cpg.coverages$V5 <- as.numeric(cpg.coverages$V5)
#colnames(cpg.coverages) <- c('CpG Island Density per mb','Proportion of Accelerated Regions displaying BGC','Species')

species_v <- c(hg19='Human',panTro4='Chimp',rheMac3='Macaque',mm10='Mouse',rn5='Rat',cavPor3='Guinea Pig',susScr3='Pig',
               bosTau7='Cow',turTru2='Dolphin',canFam3='Dog',felCat5='Cat',equCab3='Horse',myoLuc2='Brown Bat',
               loxAfr3='Elephant')

cpg.coverages$sp <-  species_v[cpg.coverages$V4]
ggplot(cpg.coverages,aes(V2,V3,label=sp)) + geom_point() + geom_smooth(method='lm', formula= y~x) + geom_text(hjust=-0.15,size=5) +
  theme_bw() + xlab('CpG Island Density per 100kb') + ylab('Proportion of Accelerated Regions displaying BGC') + theme(text = element_text(size = 20)) 
ggplot(cpg.coverages,aes(V1,V3,label=sp)) + geom_point() + geom_smooth(method='lm', formula= y~x) + geom_text(hjust=0.5,vjust=-0.75,size=5) +
  theme_bw() + xlab('Total number of CpG Islands') + ylab('Proportion of Accelerated Regions displaying BGC') + theme(text = element_text(size = 20)) 

cor(cpg.coverages$V1,cpg.coverages$V3)
cor(cpg.coverages$V1,cpg.coverages$V5)

#cpg.coverages <- cpg.coverages[cpg.coverages$V4 != 'turTru2' & cpg.coverages$V4 != 'myoLuc2' & cpg.coverages$V4 != 'felCat5' & cpg.coverages$V4 != 'cavPor3'  & cpg.coverages$V4 != 'susScr3',]
cpg.coverages$V5[cpg.coverages$V4 == 'felCat5'] <- 40
cpg.coverages$V5[cpg.coverages$V4 == 'susScr3'] <- 40
#cpg.coverages$V5[cpg.coverages$V4 == 'bosTau7'] <- 60

ggplot(cpg.coverages,aes(V5,V1,label=sp)) + geom_point() + geom_smooth(method='lm', formula= y~x) + geom_text(hjust=0.5,vjust=-0.75,size=5) +
  theme_bw() + xlab('chromosome number') + ylab('Total number of CpG Islands') + theme(text = element_text(size = 20)) 
ggplot(cpg.coverages,aes(V5,V3,label=sp)) + geom_point() + geom_smooth(method='lm', formula= y~x) + geom_text(hjust=0.5,vjust=-0.75,size=5) +
  theme_bw() + xlab('chromosome number') + ylab('Proportion of Accelerated Regions displaying BGC') + theme(text = element_text(size = 20)) 


#rat not work as don't have the gc-bias

#approx c(canFam3=25.7,bosTau7=17.8,equCab2=18.2,hg19=11.6,panTro4=11.5,mm10=9.3,rn5=9.2)
####### quickly checking whether bgc cpgs are 'younger' than other's by comparing to cat cpg islands


cat.cpgs <- import.bed('~/Downloads/cat.cpg.canFam3coords.bed')

cpg.w.bgc <- cpgIslands.gr[unique(queryHits(findOverlaps(cpgIslands.gr,cf.gcb)))]
cpg.no.bgc <- cpgIslands.gr[-unique(queryHits(findOverlaps(cpgIslands.gr,cf.gcb)))]

mcat <- matrix(c(length(cpg.no.bgc[unique(queryHits(findOverlaps(cpg.no.bgc,cat.cpgs)))]),
                 length(cpg.no.bgc[-unique(queryHits(findOverlaps(cpg.no.bgc,cat.cpgs)))]),
                 length(cpg.w.bgc[unique(queryHits(findOverlaps(cpg.w.bgc,cat.cpgs)))]),
                 length(cpg.w.bgc[-unique(queryHits(findOverlaps(cpg.w.bgc,cat.cpgs)))])
                 ),ncol = 2)
#46.71% of non-bgc CpG islands are also called as CpGs in cat, compared to 35.01% of bgc-CpGs 
fisher.test(mcat) 
#no cat or other carnivora bisulfite data, but could see if bgc dog etc are 'new' cpg islands by comparing 
#at human data
