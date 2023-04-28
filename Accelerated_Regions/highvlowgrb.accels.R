#script for looking at accelerations in high vs low turnover grb 
library(ggplot2)
library(rtracklayer)
library(data.table)
library(GenomicScores)
library(regioneR)
library(tidyverse)

setwd('~/Documents/Projects/Accelerated_Regions/')
########### prep stuff
cnes <- import.bed('CNEs.July17.bed')
accel_folder <- 'generalstuff/July Stuff/Accelerated_Regions.July2017/'
#all_lrts <- read.table(paste0(accel_folder,'CollatedLRTsunfiltered.July2017.tsv'),sep='\t',header = T) #this is crap
#all_lrts.gr <- makeGRangesFromDataFrame(all_lrts,seqnames.field = 'chr',keep.extra.columns = T)
#all_lrts.gr <- all_lrts.gr[unique(queryHits(findOverlaps(all_lrts.gr,cnes)))]
all_accels <- fread("/Users/jk2014/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/perms grb stuff/grb enrichment/accelerated_table.July2017.tsv")
all_accels <- all_accels %>% separate(coords,c("seqnames", "start"),c(':'))
all_accels <- all_accels %>% separate(start,c("start","end"),c('-'))
all_accels.gr <- makeGRangesFromDataFrame(all_accels,keep.extra.columns = T)

all_accels.gr <- all_accels.gr[unique(queryHits(findOverlaps(all_accels.gr,cnes)))]
all_accels <- as.data.frame(all_accels.gr)

cnes_ages <- import.bed('generalstuff/cne_conservation/CNEs.conserved_to.bed')
cnes_ages <- cnes_ages[unique(queryHits(findOverlaps(cnes_ages,cnes)))]

accel_counts<-rowSums(all_accels[, 6:ncol(all_accels)])
accel_counts.gr <- makeGRangesFromDataFrame(all_accels,keep.extra.columns = F)
accel_counts.gr$nAcc <- accel_counts


##

#lowturn <- makeGRangesFromDataFrame(read.table('high.low.turnover/hg19_low_turnover_grbs.bed',header = T,sep='\t'))
#highturn <- makeGRangesFromDataFrame(read.table('high.low.turnover/hg19_high_turnover_grbs.bed',header = T,sep='\t'))
lowturn <- import.bed('low_turn_hg19_canFam3_grbs.bed')
highturn <- import.bed('high_turn_hg19_canFam3_grbs.bed')
##### analysis wooo

##### Beginning- fishers tests
fishers.grbhighvlow <- function(acc.counts,n,typeof='greater'){
  if(typeof=='greater'){
    lowTurn.accs <- numOverlaps(acc.counts[acc.counts$nAcc >= n],lowturn)
    highTurn.accs <- numOverlaps(acc.counts[acc.counts$nAcc >= n],highturn)
  }
  else{
    lowTurn.accs <- numOverlaps(acc.counts[acc.counts$nAcc == n],lowturn)
    highTurn.accs <- numOverlaps(acc.counts[acc.counts$nAcc == n],highturn)
    
  }
  
  lowTurn.none <- numOverlaps(acc.counts[acc.counts$nAcc == 0],lowturn)
  highTurn.none <- numOverlaps(acc.counts[acc.counts$nAcc == 0],highturn)
  
  z <- matrix(c(lowTurn.accs,lowTurn.none,highTurn.accs,highTurn.none),ncol = 2)
  return(list(fisher.test(z),z))
}

fishers.grbhighvlow(accel_counts.gr,2)
fishers.grbhighvlow(accel_counts.gr,3)

#


accel_counts.gr$GRB <- 'None'
accel_counts.gr[unique(queryHits(findOverlaps(accel_counts.gr,lowturn)))]$GRB <- 'Low.Turnover'
accel_counts.gr[unique(queryHits(findOverlaps(accel_counts.gr,highturn)))]$GRB <- 'High.Turnover'

ac <- as.data.frame(accel_counts.gr) #want to set for each grb
###
accels.placental <- accel_counts.gr[unique(queryHits(findOverlaps(accel_counts.gr,cnes_ages[cnes_ages$score == 0])))]
accels.mammal <- accel_counts.gr[unique(queryHits(findOverlaps(accel_counts.gr,cnes_ages[cnes_ages$score == 1])))]
accels.amniote <- accel_counts.gr[unique(queryHits(findOverlaps(accel_counts.gr,cnes_ages[cnes_ages$score == 2])))]
accels.vert <- accel_counts.gr[unique(queryHits(findOverlaps(accel_counts.gr,cnes_ages[cnes_ages$score == 3])))]

fishers.grbhighvlow(accels.placental,2) #each gets stronger more accels
fishers.grbhighvlow(accels.mammal,2) #note 3 is significant here
fishers.grbhighvlow(accels.amniote,2)
fishers.grbhighvlow(accels.vert,2)

##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
#be careful here the split is working fine
plot.acc.grbs <- function(acc.counts){
  dogs <- acc.counts #[unique(queryHits(findOverlaps(acc.counts,cnes_ages[cnes_ages$score == 1])))]
  
  link.grbs.accels <- function(grbs){
    cnes.in.grb <- dogs[unique(subjectHits(findOverlaps(grbs,dogs)))]
    if(length(cnes.in.grb$nAcc) ==0){
      return(0)
    }
    return(sum(cnes.in.grb$nAcc)/length(cnes.in.grb$nAcc))
  }
  
  #link.grbs.accels <- function(grbs){
  #  cnes.in.grb <- dogs[unique(subjectHits(findOverlaps(grbs,dogs)))]
  #  return(length(cnes.in.grb[cnes.in.grb$nAcc > 0])/length(cnes.in.grb$nAcc))
  #}
  
  
  lowturn.grl <- split(lowturn, as.factor(lowturn))
  highturn.grl <- split(highturn, as.factor(lowturn))
  
  
  accels.per.grb.lowturn <- sapply(lowturn.grl, link.grbs.accels)
  accels.per.grb.highturn <- sapply(highturn.grl, link.grbs.accels)
  
  
  accels.per.grb.highturn <- data.frame(accels.per.grb.highturn,'High Turnover') ; colnames(accels.per.grb.highturn) <- c('n','Type')
  accels.per.grb.lowturn <- data.frame(accels.per.grb.lowturn,'Low Turnover') ; colnames(accels.per.grb.lowturn) <- c('n','Type')
  
  x <- rbind(accels.per.grb.highturn,accels.per.grb.lowturn)
  return(x)
}

acc.pergrb <- plot.acc.grbs(accel_counts.gr)

ggplot(acc.pergrb,aes(Type,n,fill=Type)) + geom_boxplot(outlier.alpha = 0.1) + theme_classic() + xlab('GRBs') +
  ylab('Total Acceleration Events on CNEs in GRB / Total number CNEs in GRB') + ggsci::scale_fill_jco() +
  ggpubr::stat_compare_means() + theme(text = element_text(size = 15)) 

#
acc.pl <- plot.acc.grbs(accels.placental)
acc.pl$Age <- 'Placental'
acc.mm <- plot.acc.grbs(accels.mammal)
acc.mm$Age <- 'Mammalian'
acc.an <- plot.acc.grbs(accels.amniote)
acc.an$Age <- 'Amniote'
acc.v <- plot.acc.grbs(accels.vert)
acc.v$Age <- 'Vertebrate'

age.vs.t.acc <- rbind(acc.pl,acc.mm,acc.an,acc.v)
age.vs.t.acc$Age <- factor(age.vs.t.acc$Age,levels=c('Placental','Mammalian','Amniote','Vertebrate'))
ggplot(age.vs.t.acc,aes(Age,n,fill=Type)) + geom_boxplot(outlier.alpha = 0.1) + theme_classic() + xlab('Considered CNEs conserved to..') +
  ylab('Total Acceleration Events on CNEs in GRB / Total number CNEs in GRB') +
  ggsci::scale_fill_jco()  + theme(text = element_text(size = 15)) +
  coord_cartesian(ylim=c(0,3.5))#+ ggpubr::stat_compare_means()

ggplot(age.vs.t.acc,aes(Age,n,fill=Type)) + geom_boxplot(outlier.alpha = 0.1) + theme_classic() + xlab('Considered CNEs conserved to..') +
  ylab('Total Acceleration Events on CNEs in GRB / Total number CNEs in GRB') +
  ggsci::scale_fill_jco()  + theme(text = element_text(size = 15)) + ggpubr::stat_compare_means()


#+ ylab('Proportion of CNEs with acceleration on >=1 branch')
#what if high turnover just have few placental cnes?


#

or.df <- read.table('Thesis/Figures/low.v.high.age.ors.txt',sep='\t',header = T)

or.df$Elementage <- factor(or.df$Elementage, levels=c('Placental','Mammalian','Amniote','Vertebrate'))
or.df$BranchesAccelerated <- factor(or.df$BranchesAccelerated, levels=c(1,2,3,4))

ggplot(or.df,aes(Elementage,OddsRatio,color=BranchesAccelerated)) + geom_point(size=4) + ggsci::scale_color_jco() +
  theme_bw() + theme(text = element_text(size = 15)) + xlab('CNE conserved to...') + ylab('Low vs. High Turnover GRB odds ratio')

ggplot(or.df,aes(Elementage,OddsRatio,color=BranchesAccelerated,size=P.value)) + geom_point() + ggsci::scale_color_jco() +
  theme_bw() + theme(text = element_text(size = 15)) + xlab('CNE conserved to...') + ylab('Low vs. High Turnover GRB odds ratio')

#ggplot(or.df,aes(Elementage,OddsRatio,color=BranchesAccelerated,size=P.value)) + geom_point() + scale_color_gradient(low = "red", high = "blue")  +
#  theme_bw() + theme(text = element_text(size = 15)) + xlab('CNE conserved to...') + ylab('Low vs. High Turnover GRB odds ratio')
