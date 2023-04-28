#put this in the section about functionality- asking if multi-acc regions are degrading regions
cnes <- read.table('CNEs.July17.bed')
cnes.gr <- import.bed('CNEs.July17.bed')
accel_beds <- Sys.glob('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/noGCbias/*bed')
del_beds <- Sys.glob('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/deletions/*bed')

acc <- cnes[,c(1,2,3)]
for(i in accel_beds){
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  acc[,sp] <- 0
  accel.bed <- import.bed(i)
  #
  
  acc[unique(queryHits(findOverlaps(cnes.gr,accel.bed))),][,sp] <- 1

}
dels <- cnes[,c(1,2,3)]
for(i in del_beds){
  sp <- strsplit(i,"\\/")
  sp <- sp[[1]][length(sp[[1]])]
  sp <- strsplit(sp,"\\.")[[1]][1]
  
  dels[,sp] <- 0
  del.bed <- import.bed(i)
  #
  
  dels[unique(queryHits(findOverlaps(cnes.gr,del.bed))),][,sp] <- 1
  
}

acc.sums <- cbind(dels[,c(1,2,3)],rowSums(acc[,-c(1,2,3)])) ; colnames(acc.sums) <- c('seqnames','start','end','n.Accs')
del.sums <- cbind(dels[,c(1,2,3)],rowSums(dels[,-c(1,2,3)])) ; colnames(del.sums) <- c('seqnames','start','end','n.Dels')
acc.sums.gr <- makeGRangesFromDataFrame(acc.sums,keep.extra.columns = T)
del.sums.gr <- makeGRangesFromDataFrame(del.sums,keep.extra.columns = T)

acc.sums.gr[acc.sums.gr$n.Accs >5]$n.Accs <- 5
del.sums.gr[del.sums.gr$n.Dels >5]$n.Dels <- 5

test <- data.frame(dels=del.sums.gr$n.Dels,accs=acc.sums.gr$n.Accs)
test$accs <- as.factor(test$accs)

#ggplot(test,aes(dels,fill=accs)) + geom_histogram(aes(y = ..density..)) + facet_wrap(.~accs) +
#  theme_bw() + ggsci::scale_fill_jco()

ggplot(test,aes(dels,fill=accs)) + geom_histogram(aes(y = ..density..),position='fill')  +
  theme_bw() + ggsci::scale_fill_jco() + xlab('Number of Lineages undergoing acceleration') +
  ylab('Proportion of CNEs') + guides(fill=guide_legend(title="Number of lineages with any deletion"))

test <- data.frame(dels=del.sums.gr$n.Dels,accs=acc.sums.gr$n.Accs)
test$dels <- as.factor(test$dels)

#ggplot(test,aes(accs,fill=dels)) + geom_histogram(aes(y = ..density..)) + facet_wrap(.~dels) +
#  theme_bw() + ggsci::scale_fill_jco()

ggplot(test,aes(accs,fill=dels)) + geom_histogram(aes(y = ..density..),position='fill')  +
  theme_bw() + ggsci::scale_fill_jco() + xlab('Number of Lineages with any deletion') +
  ylab('Proportion of CNEs') + guides(fill=guide_legend(title="Number of lineages undergoing acceleration"))

#get_heatmap <- function(){
#  for(a in 0:5){
#    for() d in 0:5){
#      aa <- acc.sums.gr[acc.sums.gr$n.Accs == a]
#      dd <- del.sums.gr[del.sums.gr$n.Dels == d]
#      length(unique(queryHits(findOverlaps(aa,dd))))
#    }
#  }
#}

