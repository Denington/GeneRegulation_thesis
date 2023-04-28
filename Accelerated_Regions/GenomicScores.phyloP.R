#trying out Genomic scores before i move to main script
# i am aiming on removing thius and replacing it with my own phyloP run we'll see....
#library(genomation)
library(GenomicScores)
library(data.table)
library(tidyverse)
library(magrittr)

all_accels <- fread("/Users/jk2014/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/perms grb stuff/grb enrichment/accelerated_table.July2017.tsv")
all_accels <- all_accels %>% separate(coords,c("seqnames", "start"),c(':'))
all_accels <- all_accels %>% separate(start,c("start","end"),c('-'))
all_accels.gr <- makeGRangesFromDataFrame(all_accels,keep.extra.columns = T)

accel_counts<-rowSums(all_accels[, 4:ncol(all_accels)])
accel_counts.gr <- makeGRangesFromDataFrame(all_accels,keep.extra.columns = F)
accel_counts.gr$nAcc <- accel_counts


availableGScores()

placental.pc <- getGScores("phastCons46wayPlacental.UCSC.hg19")
allsp.pc <- getGScores("phastCons100way.UCSC.hg19")

allsp.pp <- getGScores("phyloP100way.UCSC.hg19")
fitCons <- getGScores("fitCons.UCSC.hg19")

#one.acc <- gscores(placental.pc,accel_counts.gr[accel_counts.gr$nAcc == 4])

get_gscores <- function(scoretype,n){
  out <- vector(mode = "list", length = n +1)
  for(i in 0:n){
    if(i == n){
    out[i+1] <- gscores(scoretype,accel_counts.gr[accel_counts.gr$nAcc >= i])
      
    }
    else{
    out[i+1] <- gscores(scoretype,accel_counts.gr[accel_counts.gr$nAcc == i])}
  }
  return(out)
}

placental.pc.scores <- get_gscores(placental.pc,5)
ppc.df <- as.data.frame(c(placental.pc.scores[[1]],placental.pc.scores[[2]],placental.pc.scores[[3]],
                          placental.pc.scores[[4]],placental.pc.scores[[5]]))
ppc.df$nAcc <- as.character(ppc.df$nAcc)
ggplot(ppc.df,aes(nAcc,default,fill=nAcc)) + geom_boxplot(outlier.alpha = 0.01)  + theme_classic() + ggsci::scale_fill_jco() +
  ylab('Placental phastCons scores per CNE') + xlab('Number of detected Accelerated Branches on CNE') +  theme(text = element_text(size = 20)) 

vert.pc.scores <- get_gscores(allsp.pc,5)

vpc.df <- as.data.frame(c(vert.pc.scores[[1]],vert.pc.scores[[2]],vert.pc.scores[[3]],
                          vert.pc.scores[[4]],vert.pc.scores[[5]]))
vpc.df$nAcc <- as.character(ppc.df$nAcc)
ggplot(vpc.df,aes(nAcc,default,fill=nAcc)) + geom_boxplot(outlier.alpha = 0.01)  + theme_classic() + ggsci::scale_fill_jco() +
  ylab('Placental phastCons scores per CNE') + xlab('Number of detected Accelerated Branches on CNE') +  theme(text = element_text(size = 20)) 


