library(ggplot2)
library(data.table)
library(reshape2)
library(plyr)

library(GenomicRanges)
cnes <- rtracklayer::import.bed('~/Documents/Projects/Accelerated_Regions/CNEs.July17.bed') # for filtering

### modifying so can use instead of the regioneR maybe... or at least see it looks the same

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
#needs to be 5% when fewer than expected

#observed accels

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

realout<-list()
for(i in c(0, 1, 2, 3, 4, 5)){
  if(i == 5){
    realout[[as.character(i)]]<-length(real_counts[real_counts>=i])
  } else {
    realout[[as.character(i)]]<-length(real_counts[real_counts==i])
  }
}

## results of 1000 permutations of label swapping

james_shuffles <- function(accels,ntimes){
  out_list <- vector(mode = "list", length = ntimes)
  for(i in 1:ntimes){

    sim_accels <- real_values[, lapply(.SD, sample)]
    simCounts <-rowSums(sim_accels[, 2:ncol(sim_accels)])
    out_list[[i]] <- simCounts
  }
  return(out_list)
}
#jamesshuffles is added so i can look at other column combos
n.shuffs <- 1000
random_vals <- james_shuffles(real_values,n.shuffs)

#random_vals<-readRDS("/Users/jk2014/Downloads/thanksNash.JamesMultiAccels.rds")

plotData<-lapply(random_vals, function(x){
  out<-list()
  for(i in c(0, 1, 2, 3, 4, 5)){
    if(i == 5){
      out[[as.character(i)]]<-length(x[x>=i])
    } else {
      out[[as.character(i)]]<-length(x[x==i])
    }
  }
  return(do.call(cbind, out))
})

plotData<-do.call(rbind, plotData)
colnames(plotData)<-c("0", "1", "2", "3", "4", "> 4")

plotDataTest<-melt(plotData, varnames = 0)
plotDataTest[,1]<-NULL
colnames(plotDataTest)<-c("no.species", "perm.no")

plotDataTest$obs.no<-c(rep(realout[[1]], n.shuffs), rep(realout[[2]], n.shuffs), rep(realout[[3]], n.shuffs), rep(realout[[4]], n.shuffs), rep(realout[[5]], n.shuffs), rep(realout[[6]], n.shuffs))

p <- ggplot(plotDataTest, aes(x=no.species, y=perm.no)) + geom_violin() + 
  facet_wrap(~no.species, scales = "free") + theme_bw() + geom_point(aes(y=obs.no, x=no.species), colour="green") +
  stat_summary(size=5,shape=95,fun=ninety_5_percent, geom="point", colour="red") + 
  stat_summary(size=5,shape=95,fun=five_percent, geom="point", colour="red") + ylab('Observed vs. Expected number of CNEs') +
  xlab('Number of RIAMLs on CNE') + theme(text = element_text(size = 15)) 
p


p <- ggplot(plotDataTest, aes(x=perm.no)) + geom_density() + 
  facet_wrap(~no.species, scales = "free") + theme_bw() + geom_vline(aes(xintercept=obs.no), colour="green") +
  geom_vline(data=ddply(plotDataTest, "no.species", summarize, lel=ninety_5_percent(perm.no)), aes(xintercept=lel), colour="red") #+
  #geom_vline(data=ddply(plotDataTest, "no.species", summarize, lel=five_percent(perm.no)), aes(xintercept=lel), colour="red")


