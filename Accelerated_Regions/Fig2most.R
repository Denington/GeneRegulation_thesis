##### Fig 2
library(ggplot2)
library(rtracklayer)
library(data.table)
library(GenomicScores)
library(ape)
library(plyr)
library(regioneR)
library(tidyverse)
setwd('~/Documents/Projects/Accelerated_Regions/')
###########
cnes <- import.bed('CNEs.July17.bed')
accel_folder <- 'generalstuff/July Stuff/Accelerated_Regions.July2017/'
#all_lrts <- read.table(paste0(accel_folder,'CollatedLRTsunfiltered.July2017.tsv'),sep='\t',header = T) #this is crap
#all_lrts.gr <- makeGRangesFromDataFrame(all_lrts,seqnames.field = 'chr',keep.extra.columns = T)
#all_lrts.gr <- all_lrts.gr[unique(queryHits(findOverlaps(all_lrts.gr,cnes)))]
all_accels <- fread("/Users/jk2014/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/perms grb stuff/grb enrichment/accelerated_table.July2017.tsv")
all_accels <- all_accels %>% separate(coords,c("seqnames", "start"),c(':'))
all_accels <- all_accels %>% separate(start,c("start","end"),c('-'))
all_accels.gr <- makeGRangesFromDataFrame(all_accels,keep.extra.columns = T)

all_accels <- all_accels[unique(queryHits(findOverlaps(all_accels.gr,cnes.gr))),]
all_accels.gr <- all_accels.gr[unique(queryHits(findOverlaps(all_accels.gr,cnes.gr)))]
#this is because all_accels isn't fully filtered- not blatted back and forth

treemod <- read.tree('~/Downloads/Moving_stuff_to_laptop/Accelerated_Regions/scripts/hg19.100way.phyloP100way.mod')
branch.distances <- cophenetic(treemod)
#######General functions ######
ninety_5_percent<-function(x){
  out<-quantile(x, 0.95)
  names(out)<-"95%"
  return(out)
}
five_percent<-function(x){
  out<-quantile(x, 0.05)
  names(out)<-"5%"
  return(out)
}

plotPT <- function(pt){
  x <- data.frame(permuted = pt$numOverlaps$permuted,observed = pt$numOverlaps$observed)
  x$no.species <- 1
  if(mean(x$permuted) < pt$numOverlaps$observed){
  ggplot(x, aes(x=permuted)) + geom_density() + 
    theme_bw() + geom_vline(aes(xintercept=observed), colour="green") +
    geom_vline(data=ddply(x, "no.species", summarize, lel=ninety_5_percent(permuted)), aes(xintercept=lel), colour="red") #+
  }
  else{
    ggplot(x, aes(x=permuted)) + geom_density() + 
      theme_bw() + geom_vline(aes(xintercept=observed), colour="green") +
      geom_vline(data=ddply(x, "no.species", summarize, lel=five_percent(permuted)), aes(xintercept=lel), colour="red") #+
    
  }
}

lrts_to_bed <- function(elements){
  eles <- read.table(elements,header = F,sep='\t',
                col.names = c('chr','start','end','name','null_scale','alt_scale','alt_subscale','lnlratio','pval'))
  
  eles[eles$chr == 'chr1_1' | eles$chr == 'chr1_2' | eles$chr == 'chr1_3' ,]$chr <- 'chr1'
  eles[eles$chr == 'chr2_1' | eles$chr == 'chr2_2' | eles$chr == 'chr2_3' ,]$chr <- 'chr2'
  eles[eles$chr == 'chr3_1' | eles$chr == 'chr3_2' | eles$chr == 'chr3_3' ,]$chr <- 'chr3'
  
  eles.gr <- makeGRangesFromDataFrame(eles, keep.extra.columns = T, seqnames.field = 'chr',start.field = 'start',end.field = 'end')
  eles.gr <- sort(eles.gr[unique(queryHits(findOverlaps(eles.gr,cnes)))])
  return(eles.gr)
}


species_v <- c(hg19='Human',panTro4='Chimp',rheMac3='Macaque',mm10='Mouse',rn5='Rat',cavPor3='Guinea Pig',oryCun2='Rabbit',susScr3='Pig',
               bosTau7='Cow',turTru2='Dolphin',orcOrc1='Killer Whale',canFam3='Dog',felCat5='Cat',lepWed1='Seal',equCab2='Horse',myoLuc2='Brown Bat',
               pteVam1='Vampire Bat',loxAfr3='Elephant')

species_v <- c(hg19='Human',panTro4='Chimp',rheMac3='Macaque',mm10='Mouse',rn5='Rat',cavPor3='Guinea Pig',susScr3='Pig',
               bosTau7='Cow',turTru2='Dolphin',canFam3='Dog',felCat5='Cat',equCab2='Horse',myoLuc2='Brown Bat',
               loxAfr3='Elephant')

#can use this to index species_v['hg19'] etc to get actual names for plots
########## (A) Accelerated regions per species ##########  need to remove bgc from bat guinea pig and dolphin
##########  ##########  ##########  ##########  #########

branch.vs.accels <- read.table('generalstuff/July Stuff/accels.vs.branchlength.txt',sep='\t')
colnames(branch.vs.accels) <- c('Species','Branch.Length','Number.Accs')
branch.vs.accels <- branch.vs.accels #[-c(7,10,12,15,16,18,20),] 
#ggplot(branch.vs.accels,aes(Branch.Length,Number.Accs,label=Species)) + geom_point() + theme_bw() +
#  geom_text(nudge_x = 0.15) + xlab('Branch length') + ylab('Number of Accelerated Regions detected')

ggplot(branch.vs.accels,aes(Branch.Length,Number.Accs,label=Species)) + geom_smooth(method = lm) +
  geom_point() + theme_bw() +
  geom_text(hjust = 0, nudge_x = 0.0015,check_overlap = T,size=7) + xlab('Branch length') + 
  ylab('Number of Accelerated Regions detected') + theme(text = element_text(size = 15)) 

#need to make Species- num accles table using branch vs. accels  . Also need to change the species used here- why orangutan?
branch.vs.accels[,c(-2)]
cor(branch.vs.accels$Branch.Length,branch.vs.accels$Number.Accs)
########## (B) Overlaps between species ########## 
##########  ##########  ##########  ##########  #########


mutual_accels <- data.frame(matrix(0,length(species_v),length(species_v)))
colnames(mutual_accels) <- species_v
for(sp in names(species_v)){
  sp_v <- c()
  for(sp2 in names(species_v)){
    accel_in_both <- all_accels[as.vector(all_accels[,..sp] == 1) & as.vector(all_accels[,..sp2] == 1),] #fucking data.table -_-
    n.overlaps <- dim(accel_in_both)[[1]]
    sp_v <- c(sp_v,n.overlaps)
  }
  mutual_accels[,species_v[sp]] <- sp_v
}
row.names(mutual_accels) <- species_v
#now the same but with oe- can do thesame loop, but as a regioneR
mutual_accelsOE <- data.frame(matrix(0,length(species_v),length(species_v)))
colnames(mutual_accelsOE) <- species_v
mutual_accelspv <- data.frame(matrix(0,length(species_v),length(species_v)))
colnames(mutual_accelspv) <- species_v
for(sp in names(species_v)){
  sp_v <- c()
  sp_p <- c()
  
  for(sp2 in names(species_v)){
    setA <- all_accels.gr[as.vector(all_accels[,..sp] == 1)]
    setB <- all_accels.gr[as.vector(all_accels[,..sp2] == 1)]
    pt <- permTest(A= setA, B=setB, universe=all_accels.gr,ntimes = 250,randomize.function = resampleRegions,evaluate.function = numOverlaps,
                   force.parallel = T)
    oe<- pt$numOverlaps$observed / mean(pt$numOverlaps$permuted)
    sp_v <- c(sp_v,oe)
    sp_p <- c(sp_p,pt$numOverlaps$pval)
  }
  mutual_accelsOE[,species_v[sp]] <- sp_v
  mutual_accelspv[,species_v[sp]] <- sp_p
}
row.names(mutual_accelsOE) <- species_v
row.names(mutual_accelspv) <- species_v
mutual_accelspv <- as.matrix(mutual_accelspv)
diag(mutual_accelspv) <- 1

#mutual_accelsOE <- as.matrix( Matrix::forceSymmetric(mutual_accelsOE,uplo="L"))

mutual_accelsOE <- as.matrix(mutual_accelsOE)
diag(mutual_accelsOE) <- 1

col_fun = colorRamp2(c(0.5, 1, 3), c("blue", "white", "red"))
#ComplexHeatmap::Heatmap(mutual_accelsOE,cluster_rows = F,cluster_columns = F)

Heatmap(mutual_accelsOE,cluster_rows = F,cluster_columns = F, column_gap = unit(0, "mm"), border = TRUE,
        rect_gp = gpar(col = "white", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", mutual_accelsOE[i, j]), x, y, gp = gpar(fontsize = 10))},
        col = col_fun)

#yes im aware this should be vectorised. and also doing 2x too many
#also make this into a blue-red heatmap 6,7,12,13,17,18,21,22,23,25,26,27,28
mutual_accelsOE[lower.tri(mutual_accelsOE)] <- 0
mutual_accelspv[upper.tri(mutual_accelspv)] <- 0

oe.p.mutual.accels <- mutual_accelspv + mutual_accelsOE
diag(oe.p.mutual.accels) <- 1

write.table(oe.p.mutual.accels,'Thesis/Figures/oe.table.tsv',sep='\t',quote=F,row.names=T,col.names=T)
##
all_accels.sub <- all_accels #[,-c(6,7,12,13,17,18,21,22,23,25,26,27,28)]

species_pairs <- combn(species_v, 2,simplify = F)
distance.overlaps.df <- data.frame(matrix(ncol = 4, nrow = 0))  
for(i in 1:length(species_pairs)){
    evo.dist <- branch.distances[names(species_pairs[[i]][1]),names(species_pairs[[i]][2])]
    evo.row <- c(species_pairs[[i]][1][[1]] , species_pairs[[i]][2][[1]], evo.dist,mutual_accelsOE[species_pairs[[i]][1][[1]],species_pairs[[i]][2][[1]]])
    distance.overlaps.df <- rbind(distance.overlaps.df,evo.row)
}
colnames(distance.overlaps.df) <- c('Species1','Species2','Evolutionary.Distance','coacc.enrichment' )
distance.overlaps.df$Evolutionary.Distance <- as.numeric(distance.overlaps.df$Evolutionary.Distance)
distance.overlaps.df$coacc.enrichment <- as.numeric(distance.overlaps.df$coacc.enrichment)

dogg <- ggplot(distance.overlaps.df,aes(Evolutionary.Distance,coacc.enrichment)) +  geom_smooth()  + geom_point() + theme_bw() +
  xlab('Evolutionary distance between species') + ylab('Enrichment for both lineages accelerated at same elements') +
  theme(text = element_text(size = 11)) 
ggsave('Thesis/Figures/2ndpassfigures/3.4.1d.pdf',dogg,dpi = 800)

cor.test(distance.overlaps.df$Evolutionary.Distance,distance.overlaps.df$coacc.enrichment,method = 'spearman')  
  
############# plotting several examples of eg human-chimp overlaps. mod default regioneR plot
human_accels <- all_accels.gr[all_accels.gr$hg19 ==1]
chimp_accels <- all_accels.gr[all_accels.gr$panTro4 ==1]
mouse_accels <- all_accels.gr[all_accels.gr$mm10 ==1]
horse_accels <- all_accels.gr[all_accels.gr$equCab2 ==1]

plot(permTest(A=cnes[unique(queryHits(findOverlaps(cnes,chimp_accels)))],B=cnes[unique(queryHits(findOverlaps(cnes,human_accels)))],
              randomize.function = resampleRegions,evaluate.function = numOverlaps,universe=cnes,ntimes=10000))

plot(permTest(A=cnes[unique(queryHits(findOverlaps(cnes,human_accels)))],B=cnes[unique(queryHits(findOverlaps(cnes,horse_accels)))],
              randomize.function = resampleRegions,evaluate.function = numOverlaps,universe=cnes,ntimes=10000))

########## (C)  ########## 
##########  ##########  ##########  ##########  #########
#### element ages

cnes_ages <- import.bed('generalstuff/cne_conservation/CNEs.conserved_to.bed')
cnes_ages <- cnes_ages[unique(queryHits(findOverlaps(cnes_ages,cnes)))]

cnes_ages.df <- as.data.frame(cnes_ages)[,c('width','name','score')]
colnames(cnes_ages.df) <- c('width','No.Accelerated.Lineages','ConstoNum')
cnes_ages.df$Conserved.to <- 'Placental'
cnes_ages.df[cnes_ages.df$ConstoNum == 1,]$Conserved.to <- 'Mammalian'
cnes_ages.df[cnes_ages.df$ConstoNum == 2,]$Conserved.to <- 'Amniote'
cnes_ages.df[cnes_ages.df$ConstoNum == 3,]$Conserved.to <- 'Vertebrate'
cnes_ages.df$ConstoNum <- NULL 

cnes_ages.df$No.Accelerated.Lineages <- as.numeric(cnes_ages.df$No.Accelerated.Lineages)
cnes_ages.df[cnes_ages.df$No.Accelerated.Lineages >= 6,]$No.Accelerated.Lineages <- 6
cnes_ages.df[cnes_ages.df$No.Accelerated.Lineages == 6,]$No.Accelerated.Lineages <- '6+'

cnes_ages.df$No.Accelerated.Lineages <- factor(cnes_ages.df$No.Accelerated.Lineages,levels= c('0','1','2','3','4','5','6+'))


cne.age.vs.accels <- cnes_ages.df %>% count(Conserved.to,No.Accelerated.Lineages)
cne.age.vs.accels$Conserved.to <- factor(cne.age.vs.accels$Conserved.to,levels = c('Placental','Mammalian','Amniote','Vertebrate'))


ggplot(cne.age.vs.accels,aes(x=No.Accelerated.Lineages,y=n,fill=Conserved.to)) + geom_bar(stat='identity',position = 'fill') +
  theme_classic() + ggsci::scale_fill_jco() + xlab('Number of Lineages Accelerated in CNE') + 
  ylab('Proportion of Tested Elements') + theme(text = element_text(size = 20)) 

ggplot(cnes_ages.df,aes(No.Accelerated.Lineages,width,fill=No.Accelerated.Lineages)) + geom_boxplot(outlier.alpha = 0.005)  +
  theme_classic() + ggsci::scale_fill_jco() + coord_cartesian(ylim=c(0,1000)) + xlab('Number of Accelerated Lineages on CNE') +
  ylab('CNE width (bp)') + theme(text = element_text(size = 20)) + ggpubr::stat_compare_means(comparisons = )
  

################### Ancestral branch - closer LRT scores?
###################
library(ggpubr)

mouse_accels <- import.bed(paste0(accel_folder,'mm10.AcceleratedRegions.bed'))
rat_accels <- import.bed(paste0(accel_folder,'rn5.AcceleratedRegions.bed'))
dog_accels <- import.bed(paste0(accel_folder,'canFam3.AcceleratedRegions.bed'))

human_accels <- import.bed(paste0(accel_folder,'hg19.AcceleratedRegions.bed'))
chimp_accels <- import.bed(paste0(accel_folder,'panTro4.AcceleratedRegions.bed'))

human_chimp_anc_accels <- import.bed(paste0(accel_folder,'hg19-panTro4.AcceleratedRegions.bed'))

#hum_chimpacels <- human_accels[unique(queryHits(findOverlaps(human_accels,chimp_accels)))]
#hum_ancacels <- human_accels[unique(queryHits(findOverlaps(human_accels,human_chimp_anc_accels)))]
#chimp_ancacels <- chimp_accels[unique(queryHits(findOverlaps(chimp_accels,human_chimp_anc_accels)))]
#all_ancacels <- chimp_ancacels[unique(queryHits(findOverlaps(chimp_ancacels,hum_ancacels)))]


human_lrts <-  lrts_to_bed('generalstuff/LRT_files/July2017/hg19/hg19.LRT.bed')
chimp_lrts <-  lrts_to_bed('generalstuff/LRT_files/July2017/panTro4/panTro4.LRT.bed')
hc_anc_lrts <-  lrts_to_bed('generalstuff/LRT_files/July2017/hg19-panTro4/hg19-panTro4.LRT.bed')
mouse_lrts <-  lrts_to_bed('generalstuff/LRT_files/July2017/mm10/mm10.LRT.bed')
rat_lrts <-  lrts_to_bed('generalstuff/LRT_files/July2017/rn5/rn5.LRT.bed')
mr_anc_lrts <-  lrts_to_bed('generalstuff/LRT_files/July2017/mm10-rn5/mm10-rn5.LRT.bed')
dog_lrts <-  lrts_to_bed('generalstuff/LRT_files/July2017/canFam3/canFam3.LRT.bed')


human_lrts$humanLNL <- human_lrts$lnlratio
human_lrts$chimpLNL <- chimp_lrts$lnlratio
human_lrts$ancLNL <- hc_anc_lrts$lnlratio
human_lrts$mouseLNL <- mouse_lrts$lnlratio
human_lrts$ratLNL <- rat_lrts$lnlratio
human_lrts$mouseratLNL <- mr_anc_lrts$lnlratio
human_lrts$dogLNL <- dog_lrts$lnlratio


humACCs <- human_lrts[unique(queryHits(findOverlaps(human_lrts,human_accels)))]; humACCs$Accelerated <- 'Human'
chimpACCs <- human_lrts[unique(queryHits(findOverlaps(human_lrts,chimp_accels)))]; chimpACCs$Accelerated <- 'Chimp'
mouseACCs <- human_lrts[unique(queryHits(findOverlaps(human_lrts,mouse_accels)))]; mouseACCs$Accelerated <- 'Mouse'
ratACCs <- human_lrts[unique(queryHits(findOverlaps(human_lrts,rat_accels)))]; ratACCs$Accelerated <- 'Rat'
dogACCs <- human_lrts[unique(queryHits(findOverlaps(human_lrts,dog_accels)))]; dogACCs$Accelerated <- 'Dog'
ancACCs <- human_lrts[unique(queryHits(findOverlaps(human_lrts,human_chimp_anc_accels)))]; ancACCs$Accelerated <- 'Human-Chimp Ancestor'
allACCs <- c(humACCs,chimpACCs,mouseACCs,ratACCs,ancACCs,dogACCs)
noACCs <- human_lrts[-unique(queryHits(findOverlaps(human_lrts,allACCs)))] ; noACCs$Accelerated <- 'None'
accels.vs.lnl <- c(allACCs,noACCs)
#human_lrts$Accelerated <- 'None'

#human_lrts[unique(queryHits(findOverlaps(human_lrts,human_accels)))]$Accelerated <- 'Human'#human_lrts[unique(queryHits(findOverlaps(human_lrts,chimp_accels)))]$Accelerated <- 'Chimp'#human_lrts[unique(queryHits(findOverlaps(human_lrts,human_chimp_anc_accels)))]$Accelerated <- 'Ancestor'#human_lrts[unique(queryHits(findOverlaps(human_lrts,mouse_accels)))]$Accelerated <- 'Mouse'
#human_lrts[unique(queryHits(findOverlaps(human_lrts,queryHits(findOverlaps(chimp_accels,human_accels)))))]$Accelerated <- 'Human&Chimp'

hl <- as.data.frame(accels.vs.lnl)
hl$Accelerated <- factor(hl$Accelerated, levels=c('Human','Chimp','Human-Chimp Ancestor','Mouse','Rat','Dog','None'))
#ggplot(hl,aes(Accelerated,ancLNL)) + geom_violin() + coord_cartesian(ylim=c(0,5)) + theme_classic()
ggplot(hl[hl$ancLNL != 0 & hl$Accelerated != 'Human-Chimp Ancestor',],aes(Accelerated,ancLNL,fill=Accelerated)) + 
  geom_boxplot(outlier.alpha = 0.0, show.legend=F) +
  coord_cartesian(ylim=c(0,5)) + theme_classic() + ggsci::scale_fill_jco() + xlab('CNE accelerated in...') +
  ylab('Acceleration Likelihood score Human-Chimp ancestor') + theme(text = element_text(size = 10)) +
  stat_compare_means(comparisons = list(c('Chimp','Mouse'),c('Chimp','Rat'),c('Human','Mouse'),c('Human','Rat')))
#ggplot(hl[hl$humanLNL != 0 & hl$Accelerated != 'Human',],aes(Accelerated,humanLNL,fill=Accelerated)) + geom_boxplot(outlier.alpha = 0.05) + 
#  coord_cartesian(ylim=c(0,5)) + theme_classic() + ggsci::scale_fill_jco() + xlab('CNE accelerated in...') + ylab('Acceleration Likelihood score Human')
#ggplot(hl[hl$chimpLNL != 0 & hl$Accelerated != 'Chimp',],aes(Accelerated,chimpLNL,fill=Accelerated)) + geom_boxplot(outlier.alpha = 0.05) +
#  coord_cartesian(ylim=c(0,5)) + theme_classic() + ggsci::scale_fill_jco() + xlab('CNE accelerated in...') + ylab('Acceleration Likelihood score Chimp')
#ggplot(hl[hl$ancLNL != 0,],aes(Accelerated2,ancLNL)) + geom_violin() + coord_cartesian(ylim=c(0,5)) + theme_classic()

ggplot(hl[hl$mouseratLNL != 0 & hl$Accelerated != 'Human-Chimp Ancestor',],aes(Accelerated,mouseratLNL,fill=Accelerated)) + 
  geom_boxplot(outlier.alpha = 0.05) + coord_cartesian(ylim=c(0,10)) + theme_classic() + ggsci::scale_fill_jco() + 
  xlab('CNE accelerated in...') + ylab('Acceleration Likelihood score Mouse-Rat ancestor') + theme(text = element_text(size = 20)) +
  stat_compare_means(comparisons = list(c('Chimp','Mouse'),c('Chimp','Rat'),c('Human','Mouse'),c('Human','Rat')))

#ggplot(hl[hl$ratLNL != 0 & hl$Accelerated != 'Rat',],aes(Accelerated,ratLNL,fill=Accelerated)) + geom_boxplot(outlier.alpha = 0.05) + 
#  coord_cartesian(ylim=c(0,10)) + theme_classic() + ggsci::scale_fill_jco() + xlab('CNE accelerated in...') + ylab('Acceleration Likelihood score Rat')

hh <- hl[hl$humanLNL >0.1 & hl$ancLNL >0.1,]
smoothScatter(hh$ancLNL,hh$humanLNL,ylim = c(0,7.5),xlim = c(0,7.5),
      ylab = 'Likelihood score CNE is accelerated on Human Branch',xlab='Likelihood score CNE is accelerated on Human-Chimp Ancestor Branch')

#hh <- hl[hl$humanLNL >0.1 & hl$chimpLNL >0.1,]
#smoothScatter(hh$humanLNL,hh$chimpLNL,ylim = c(0,7.5),xlim = c(0,7.5),
#      xlab = 'Likelihood score CNE is accelerated on Human Branch',ylab='Likelihood score CNE is accelerated on Chimp Branch')


hh <- hl[hl$chimpLNL >0.1 & hl$ancLNL >0.1,]
smoothScatter(hh$ancLNL,hh$chimpLNL, ylim = c(0,7.5),xlim = c(0,7.5),
          ylab = 'Likelihood score CNE is accelerated on Chimp Branch',xlab='Likelihood score CNE is accelerated on Human-Chimp Ancestor Branch')


#hh <- hl[hl$humanLNL >0.1 & hl$mouseLNL >0.1,]
#smoothScatter(hh$humanLNL,hh$mouseLNL,ylim = c(0,7.5),xlim = c(0,7.5),
#              xlab = 'Likelihood score CNE is accelerated on Human Branch',ylab='Likelihood score CNE is accelerated on Mouse Branch')

#abline(a=-0, b=1,col='red') 

#
#hh <- hl[hl$mouseLNL >0.1 & hl$ratLNL >0.1,]
#smoothScatter(hh$mouseLNL,hh$ratLNL,ylim = c(0,7.5),xlim = c(0,7.5),
#      xlab = 'Likelihood score CNE is accelerated on Mouse Branch',ylab='Likelihood score CNE is accelerated on Rat Branch')

#abline(a=-0, b=1,col='red') 

hh <- hl[hl$mouseLNL >0.1 & hl$mouseratLNL >0.1,]
smoothScatter(hh$mouseratLNL,hh$mouseLNL,ylim = c(0,7.5),xlim = c(0,7.5),
    ylab = 'Likelihood score CNE is accelerated on Mouse Branch',xlab='Likelihood score CNE is accelerated on Mouse-Rat Ancestor Branch')

#abline(a=-0, b=1,col='red') 

hh <- hl[hl$ratLNL >0.1 & hl$mouseratLNL >0.1,]
smoothScatter(hh$mouseratLNL,hh$ratLNL,ylim = c(0,7.5),xlim = c(0,7.5),
        ylab = 'Likelihood score CNE is accelerated on Rat Branch',xlab='Likelihood score CNE is accelerated on Mouse-Rat Ancestor Branch')

#abline(a=-0, b=1,col='red') 

########tsars overlapping accels ######
#########################################

tsars <- unique(import.bed('generalstuff/tsar_beds/tsars.bed'))


accel_counts<-rowSums(all_accels[, 4:ncol(all_accels)])
accel_counts.gr <- makeGRangesFromDataFrame(all_accels,keep.extra.columns = F)
accel_counts.gr$nAcc <- accel_counts

tsar.0.pt <- permTest(A= accel_counts.gr[accel_counts.gr$nAcc == 0 ], B=tsars, universe=accel_counts.gr,ntimes = 1000,randomize.function = resampleRegions,evaluate.function = numOverlaps,
              force.parallel = T)
tsar.1.pt <- permTest(A= accel_counts.gr[accel_counts.gr$nAcc ==1], B=tsars, universe=accel_counts.gr,ntimes = 1000,randomize.function = resampleRegions,evaluate.function = numOverlaps,
              force.parallel = T)
tsar.2.pt <- permTest(A= accel_counts.gr[accel_counts.gr$nAcc ==2], B=tsars, universe=accel_counts.gr,ntimes = 1000,randomize.function = resampleRegions,evaluate.function = numOverlaps,
                      force.parallel = T)
tsar.3.pt <- permTest(A= accel_counts.gr[accel_counts.gr$nAcc ==3], B=tsars, universe=accel_counts.gr,ntimes = 1000,randomize.function = resampleRegions,evaluate.function = numOverlaps,
                      force.parallel = T)
tsar.4p.pt <- permTest(A= accel_counts.gr[accel_counts.gr$nAcc >=4], B=tsars, universe=accel_counts.gr,ntimes = 1000,randomize.function = resampleRegions,evaluate.function = numOverlaps,
                      force.parallel = T)

tsar.df <- data.frame( Overlaps=c(tsar.0.pt$numOverlaps$observed, mean(tsar.0.pt$numOverlaps$permuted),
                       tsar.1.pt$numOverlaps$observed, mean(tsar.1.pt$numOverlaps$permuted),
                       tsar.2.pt$numOverlaps$observed, mean(tsar.2.pt$numOverlaps$permuted),
                       tsar.3.pt$numOverlaps$observed, mean(tsar.3.pt$numOverlaps$permuted),
                       tsar.4p.pt$numOverlaps$observed, mean(tsar.4p.pt$numOverlaps$permuted)) ,
            Accels=c('0','0','1','1','2','2','3','3','4+','4+'),
            obsexp=c(rep(c('Observed','Expected'),5)))

ggplot(tsar.df,aes(Accels,Overlaps,fill=obsexp)) + geom_bar(stat='identity',position='dodge') + theme_classic() +
  ggsci::scale_fill_jco() + xlab('Number of lineages Accelerated') + 
  ylab('Accelerated Regions also accelerated in the ancestral mammal') +
  theme(text = element_text(size = 15)) 
#just put the oe/text over the plot



####################################################################################################################################
################## Just getting barcharts (again...) with the proportion of elements in eg introns #################################

introns <- import.bed('genome.features/ucsc.introns.bed')
tpUTR <- import.bed('genome.features/3pUTRs.bed')
fpUTR <- import.bed('genome.features/5pUTRs.bed')
linc <- import.bed('genome.features/lincrnas.bed')
retro <- import.bed('genome.features/retroposed.genes.bed')
snomrna <- import.bed('genome.features/snomrna.bed')
cpgs <- import.bed('genome.features/hg19.CpGislands.bed')
proms <- import.bed('genome.features/epd.new.promoters.bed')
marklist <- list('intron'=introns,'3\'UTR'=tpUTR,'5\'UTR'=fpUTR,'linc'=linc,'CpG island'=cpgs,'promoter'=proms) #,snomrna=snomrna)
#add cpgislands and promoters
make.another.barplot <- function(elements,marks,elety){
  out <- df <- data.frame(Mark=character(),
                                     Proportion=numeric(), 
                                     Number.Accelerated.Lineages=character())
  for(m in names(marks)){
  a <- length(unique(queryHits(findOverlaps(elements,marks[m][[1]]))))
  prop <- a/length(elements)
  o <- data.frame(m,prop,elety); colnames(o) <- c('Mark','Proportion','RAIML.Number')
  out <- rbind(out,o)
  }
  ig <- 1- sum(out$Proportion)
  out <- rbind(out,data.frame(Mark='Intergenic',Proportion=ig,RAIML.Number=elety))
  return(out)
}

b0 <- make.another.barplot(accel_counts.gr[accel_counts.gr$nAcc ==0],marklist,'0')
b1 <- make.another.barplot(accel_counts.gr[accel_counts.gr$nAcc ==1],marklist,'1')
b2 <- make.another.barplot(accel_counts.gr[accel_counts.gr$nAcc ==2],marklist,'2')
b3 <- make.another.barplot(accel_counts.gr[accel_counts.gr$nAcc ==3],marklist,'3')
b4 <- make.another.barplot(accel_counts.gr[accel_counts.gr$nAcc >=4],marklist,'4+')
#b5 <- make.another.barplot(accel_counts.gr[accel_counts.gr$nAcc >5],marklist,'5+')
bb <- rbind(b0,b1,b2,b3,b4)#,b5)

ggplot(bb,aes(Mark,Proportion,fill=RAIML.Number)) + geom_bar(stat='identity',position='dodge') +
  theme_classic() + ggsci::scale_fill_jco() + xlab('Feature') + theme(text = element_text(size = 20)) 

######
deletions <- read.table('generalstuff/July Stuff/deletions/mammal_deletions.tsv',sep='\t',col.names = c('seqnames','start','end','species',                                                                                                       'Accel.pval','x1','x2','x3'))
deletions.gr <- makeGRangesFromDataFrame(deletions,keep.extra.columns = T)



