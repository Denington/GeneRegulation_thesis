#####
##### for each grb, getting how many accels for each species in that grb

library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(tidyverse)
setwd('~/Documents/Projects/Accelerated_Regions/')

#lets do this as a for list of species accel files, have a set up grb 3 column df and cbind the sapply length overlaps 
cnes <- import.bed('CNEs.July17.bed')

#grbs.gr and grl are not ordered the same
#do i want to use a broader grb set? e.g  hg19.canFam3.98.50.processed.threshold.grbs.bed
x <- read.table('hg19.canFam3.98.50.processed.threshold.grbs.bed',col.names = c('seqnames','start','end','loci','score','x','start2','end2','color'),skip = 1)
#x <- read.table('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/grb_boundaries_genes.hg19.txt',col.names = c('seqnames','start','end','target_gene'))
grbs.gr <- unique(makeGRangesFromDataFrame(x,keep.extra.columns = T))
grbs.gr <- grbs.gr[width(grbs.gr) > 100000]
#grbs.gr <- grbs.gr[unique(grbs.gr$target_gene)]
grbs.gr <- grbs.gr[order(as.factor(grbs.gr$loci))]
#grbs.gr <- import.bed('hg19_GRBS.bed')
grbs.grl <- split(grbs.gr, as.factor(grbs.gr$loci))
grb <- as.data.frame(grbs.gr) #read.table('hg19_GRBS.bed',sep='\t',header = F) #want one with names. could just take 'accels_per_grb.July2017.tsv'
accel_beds <- Sys.glob('generalstuff/July Stuff/Accelerated_Regions.July2017/*bed')
for(i in accel_beds){
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  accel.bed <- import.bed(i)
  #
  link.accels.grbs <- function(grbs.grl){
    accs.in.grb <- accel.bed[unique(subjectHits(findOverlaps(grbs.grl,accel.bed)))]
    
    return(length(accs.in.grb))
  }
  #
  o <- sapply(grbs.grl,link.accels.grbs)
  grb[,sp] <- o
}

accel.bed <- cnes
grb$cnes <- sapply(grbs.grl,link.accels.grbs)


x <- gather(grb[,-c(1,2,3,4,5,6)],species,count)
#sp_keep <- c('hg19','panTro4','mm10','rn5','cavPor3','bosTau7','canFam3','equCab2')
sp_keep <- c('hg19','panTro4','rheMac3','mm10','rn5','cavPor3','bosTau7','turTru2','susScr3','equCab2','felCat5','canFam3','myoLuc2','loxAfr3')
x <- x[x$species %in% sp_keep,]
x$species <- factor(x$species,levels = sp_keep)
x$count <- as.numeric(x$count)
x$sp.names <- species_v[x$species]
x$sp.names <- factor(x$sp.names,levels=species_v)
ggplot(x[x$species != 'cnes',],aes(sp.names,count)) + geom_boxplot() + theme_classic() + 
  xlab('Species') + ylab('Number of Accelerated Regions in GRB')

#alright, that's pure number. I also want an
#for soms, try with scaling, with /total cnes and with /cnes and then scale(). Also try without human/chimp and other very low accel species

#and accel enrichment in cnes
# can i do something like grb$enrich <- fishertest(grb$sp,grb$cnes, allgrbs$sp - grb$sp, grb$cnes, allgrbs$cnes - grb$cnes)
# just re-do the above for loop, but with the link.accels.grbs instead returning a fisher test. need the total cnes/ accels in grbs for
# each species available

grb.fish.enrich <- as.data.frame(grbs.gr) #read.table('hg19_GRBS.bed',sep='\t',header = F) #want one with names. could just take 'accels_per_grb.July2017.tsv'
for(i in accel_beds){
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  accel.bed <- import.bed(i)
  #
  accels.grbs.enrich <- function(grbs.grl){
    accs.grb <- length(accel.bed[unique(subjectHits(findOverlaps(grbs.grl,accel.bed)))])
    cnes.grb <- length(cnes[unique(subjectHits(findOverlaps(grbs.grl,cnes)))])
    
    acc.allgrbs <- length(accel.bed[unique(subjectHits(findOverlaps(grbs.gr,accel.bed)))])
    cnes.allgrb <- length(cnes[unique(subjectHits(findOverlaps(grbs.gr,cnes)))])
    c <- acc.allgrbs - accs.grb ; d <- cnes.allgrb - cnes.grb
    m <- matrix(c(accs.grb,cnes.grb,c,d),ncol = 2)
    #i think the order of the matrix needs to be reversed?
    return(fisher.test(m)$estimate)  #p.value
  }
  #
  o <- sapply(grbs.grl,accels.grbs.enrich)
  grb.fish.enrich[,sp] <- o
}


y <- gather(grb.fish.enrich[,-c(1,2,3,4,5,6)],species,oe)
sp_keep <- c('hg19','panTro4','mm10','rn5','cavPor3','bosTau7','canFam3','equCab2')
y <- y[y$species %in% sp_keep,]
y$species <- factor(y$species,levels = sp_keep)
y$oe <- as.numeric(y$oe)
ggplot(y[y$species != 'cnes',],aes(species,oe)) + geom_boxplot() + theme_classic()
#
grb.fish.enrichP <- as.data.frame(grbs.gr)
for(i in accel_beds){
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  accel.bed <- import.bed(i)
  #
  accels.grbs.enrich <- function(grbs.grl){
    accs.grb <- length(accel.bed[unique(subjectHits(findOverlaps(grbs.grl,accel.bed)))])
    cnes.grb <- length(cnes[unique(subjectHits(findOverlaps(grbs.grl,cnes)))])
    
    acc.allgrbs <- length(accel.bed[unique(subjectHits(findOverlaps(grbs.gr,accel.bed)))])
    cnes.allgrb <- length(cnes[unique(subjectHits(findOverlaps(grbs.gr,cnes)))])
    c <- acc.allgrbs - accs.grb ; d <- cnes.allgrb - cnes.grb
    m <- matrix(c(accs.grb,cnes.grb,c,d),ncol = 2)
    #i think the order of the matrix needs to be reversed?
    return(fisher.test(m)$p.value)  #p.value
  }
  #
  o <- sapply(grbs.grl,accels.grbs.enrich)
  grb.fish.enrichP[,sp] <- o
}

bh.correct <- function(x){
  return(p.adjust(x,method = 'BH'))
}

grb.enrich.pfdr <-as.data.frame(apply(grb.fish.enrichP[,-c(1:11)],2,bh.correct))
grb.fish.enrichP2 <- cbind(grb.fish.enrichP[,c(1:11)],grb.enrich.pfdr)
minp <- function(x){return(min(x) < 0.05)}

enriched <- apply(grb.fish.enrichP2[,-c(1:11)],1,minp)
#enriched <- apply(grb.fish.enrichP2[,-c(1:11,17,23,26,29,33,36)],1,minp)
length(enriched[enriched==T])


accel.bed <- cnes
grb.fish.enrichP$n.cnes <- sapply(grbs.grl,link.accels.grbs)


#want species arranged here nicely
library(DescTools)
library(som)
library(kohonen)

wons <- function(x){
  return(Winsorize(x,probs = c(0,0.95)))
}
gfe <- grb.fish.enrich[,-c(1:11,17,19,23,24,25,26,28,29,30,33,35,36)]
gfe <- gfe[grb$cnes > 60,]
gfe <- apply(gfe,2,wons) #i think wins before scale but idk
gfe <- scale(gfe)

set.seed(123456)

test_som<-kohonen::som(as.matrix(gfe), rlen=1000, grid=somgrid(5,5,"hexagonal"))
#test_som<-kohonen::som(as.matrix(gfe2[,c(7,16,1,2,10,23,26)]), rlen=1000, grid=somgrid(4,5,"hexagonal"))
#test_som<-kohonen::som(as.matrix(gfe2[,c(7,16,10,21,23,1,26,4,2,5)]), rlen=1000, grid=somgrid(4,5,"hexagonal"))
#test_som<-kohonen::som(as.matrix(gfe2[,c(5,8,9,6,1,10,3,2,4,7)]), rlen=1000, grid=somgrid(4,5,"hexagonal"))

plot(test_som, palette.name=rainbow)

plot(test_som, type='counts')
plot(test_som,type="mapping")
test_som$codes
grb.fish.enrich[grb$cnes > 60,][test_som$unit.classif ==1,]$target_gene

#################
species_v <- c(hg19='Human',panTro4='Chimp',rheMac3='Macaque',mm10='Mouse',rn5='Rat',cavPor3='Guinea Pig',oryCun2='Rabbit',susScr3='Pig',
               bosTau7='Cow',turTru2='Dolphin',orcOrc1='Killer Whale',canFam3='Dog',felCat5='Cat',lepWed1='Seal',equCab2='Horse',myoLuc2='Brown Bat',
               pteVam1='Vampire Bat',loxAfr3='Elephant')

####
corvec <- c()
for(sp in colnames(grb.fish.enrich)[12:length(colnames(grb.fish.enrich))]){
  for(sp2 in colnames(grb.fish.enrich)[12:length(colnames(grb.fish.enrich))]){
    if(sp == sp2){next    }
    correl <- cor(grb.fish.enrich[,sp],grb.fish.enrich[,sp2])
    corx <- setNames(correl, paste0(species_v[sp],'.', species_v[sp2]) )
    corvec <- c(corvec,corx)
  }
}
enrich.cor <- data.frame('species' = names(corvec),'correl'=as.numeric(corvec)) 
ggplot(enrich.cor,aes(species,correl,label=species)) + geom_point() + geom_label() #how to make 1d?
######################################################################################################################################
###################################Are ARs found in clusters to the same proportion as CNEs? #########################################

all.accels.in.grbs <- all_accels.gr[rowSums(as.data.frame(all_accels.gr)[,-c(1,2,3,4,5)]) > 0,] 
aaig <- length(all.accels.in.grbs[unique(queryHits(findOverlaps(all.accels.in.grbs,grbs.gr)))]) 
naig <- length(all.accels.in.grbs[-unique(queryHits(findOverlaps(all.accels.in.grbs,grbs.gr)))]) 

no.accel.cnes.in.grbs <- all_accels.gr[rowSums(as.data.frame(all_accels.gr)[,-c(1,2,3,4,5)]) == 0,] 
ncig <- length(no.accel.cnes.in.grbs[unique(queryHits(findOverlaps(no.accel.cnes.in.grbs,grbs.gr)))]) 
ncog <- length(no.accel.cnes.in.grbs[-unique(queryHits(findOverlaps(no.accel.cnes.in.grbs,grbs.gr)))]) 

fisher.test(matrix(c(aaig,naig,ncig,ncog),ncol=2))
######################################################################################################################################
###################################DO ARs clusters to the same proportion as CNEs? #########################################

library(regioneR)

meanDistance2 <- function(A, ...) {
  
  if(!hasArg(A)) stop("A is missing")

  A <- toGRanges(A)

  d <- GenomicRanges::distanceToNearest(A)
  
  return(mean(as.matrix(d@elementMetadata@listData$distance)[,1], na.rm=TRUE)) #--> BioC 2.13 
  #return(mean(d@listData$distance, na.rm=TRUE)) #--> BioC 2.11
  
}

ptAll <- permTest(A= all_accels.gr[rowSums(as.data.frame(all_accels.gr)[,-c(1,2,3,4,5)]) > 0,], universe=all_accels.gr,ntimes = 1000,randomize.function = resampleRegions,evaluate.function = meanDistance2,
               force.parallel = T)

ptAll.mult <- permTest(A= all_accels.gr[rowSums(as.data.frame(all_accels.gr)[,-c(1,2,3,4,5)]) > 2,], universe=all_accels.gr,ntimes = 1000,randomize.function = resampleRegions,evaluate.function = meanDistance2,
               force.parallel = T)

ptMouse <- permTest(A= all_accels.gr[all_accels.gr$mm10 == 1], universe=all_accels.gr,ntimes = 1000,randomize.function = resampleRegions,evaluate.function = meanDistance2,
               force.parallel = T)

ptHum <- permTest(A= all_accels.gr[all_accels.gr$hg19 == 1], universe=all_accels.gr,ntimes = 1000,randomize.function = resampleRegions,evaluate.function = meanDistance2,
               force.parallel = T)


rbind(data.frame(distance=pt$meanDistance2$permuted,typeof='Permuted'),
data.frame(distance=as.data.frame(distanceToNearest(all_accels.gr[rowSums(as.data.frame(all_accels.gr)[,-c(1,2,3,4,5)]) > 0,]))$distance,typeof='Observed'))
#looks ugly as comparing a mean to an actual?
#ARs are closer together on average than CNEs




multAcc.gr <- makeGRangesFromDataFrame(all_accels)
multAcc.gr$nAcc <- rowSums(all_accels[,-c(1,2,3)])


accels.grbs.enrich.m <- function(grbs.grl){
    accs.grb <- length(multAcc.gr[multAcc.gr$nAcc > 3][unique(subjectHits(findOverlaps(grbs.grl,multAcc.gr[multAcc.gr$nAcc > 3])))])
    cnes.grb <- length(cnes[unique(subjectHits(findOverlaps(grbs.grl,cnes)))])
    
    acc.allgrbs <- length(multAcc.gr[multAcc.gr$nAcc > 3][unique(subjectHits(findOverlaps(grbs.gr,multAcc.gr[multAcc.gr$nAcc > 3])))])
    cnes.allgrb <- length(cnes[unique(subjectHits(findOverlaps(grbs.gr,cnes)))])
    c <- acc.allgrbs - accs.grb ; d <- cnes.allgrb - cnes.grb
    m <- matrix(c(accs.grb,cnes.grb,c,d),ncol = 2)
    #i think the order of the matrix needs to be reversed?
    return(fisher.test(m)$estimate)  #p.value
  }
  #
o <- sapply(grbs.grl,accels.grbs.enrich.m)

accels.grbs.enrich.m <- function(grbs.grl){
  accs.grb <- length(multAcc.gr[multAcc.gr$nAcc > 3][unique(subjectHits(findOverlaps(grbs.grl,multAcc.gr[multAcc.gr$nAcc > 3])))])
  cnes.grb <- length(cnes[unique(subjectHits(findOverlaps(grbs.grl,cnes)))])
  
  acc.allgrbs <- length(multAcc.gr[multAcc.gr$nAcc > 3][unique(subjectHits(findOverlaps(grbs.gr,multAcc.gr[multAcc.gr$nAcc > 3])))])
  cnes.allgrb <- length(cnes[unique(subjectHits(findOverlaps(grbs.gr,cnes)))])
  c <- acc.allgrbs - accs.grb ; d <- cnes.allgrb - cnes.grb
  m <- matrix(c(accs.grb,cnes.grb,c,d),ncol = 2)
  #i think the order of the matrix needs to be reversed?
  return(fisher.test(m)$p.value)  #p.value
}
#
op <- sapply(grbs.grl,accels.grbs.enrich.m)
op2 <- p.adjust(op,method = 'BH')
#grb.fish.enrich[,sp] <- o

#save.image(file="scripts/all.accels.linkedtoGRBS.RData") 