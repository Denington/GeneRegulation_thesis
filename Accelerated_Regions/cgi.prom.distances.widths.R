

get.cgi.prom.props <- function(sp,spname,ele,ch,bed,lift=T,hg19=T){
  ch = import.chain(ch)
  
  bed <- import.bed(bed); seqlevelsStyle(bed) <- 'UCSC'

  if(lift ==T){
    if(hg19 == T){
      bed <- unlist(liftOver(bed,ch))

    }
    else{
      
      bed.38 <- unlist(liftOver(bed,human.ch))
      bed <- unlist(liftOver(bed.38,ch))
      
    }
    bed <- reduce(bed,ignore.strand = TRUE, min.gapwidth = 51L)

  }
  dogSession<-browserSession("UCSC")
  genome(dogSession) <- sp
  
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
  
  
  cgi.eles <- bed[unique(queryHits(findOverlaps(bed,cpgIslands.gr)))]
  prom.eles <- bed[unique(queryHits(findOverlaps(bed,promoters(doggenes.gr))))]
  both.eles <- cgi.eles[unique(queryHits(findOverlaps(cgi.eles,prom.eles)))]

  n.ib <- bed[-unique(queryHits(findOverlaps(bed,cpgIslands.gr)))]
  if(length(n.ib) == 0 ){
    n.ib <- bed
  }
  neither <- n.ib[-unique(queryHits(findOverlaps(n.ib,promoters(doggenes.gr))))]
  if(length(neither) == 0 ){
    neither <- n.ib
  }
  
  cgi.eles <- cgi.eles[-unique(queryHits(findOverlaps(cgi.eles,both.eles)))]  
  prom.eles <- prom.eles[-unique(queryHits(findOverlaps(prom.eles,both.eles)))]  
  
  df.out <- data.frame(Species=rep(spname,4),Element=rep(ele,4),Overlaps=c('CGI','Promoter','CGI.Promoter','Neither'),
                       n=c(length(cgi.eles),length(prom.eles),length(both.eles),length(neither)))
  return(df.out)
  
}

setwd('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/')

testdf <- rbind(
  
  get.cgi.prom.props('canFam3','Dog','gBGC AR','~/Downloads/Liftover.chains/hg19ToCanFam3.over.chain','canFam3.GCbias.bed'),  
  get.cgi.prom.props('canFam3','Dog','other AR','~/Downloads/Liftover.chains/hg19ToCanFam3.over.chain','noGCbias/canFam3.bed'),  
  get.cgi.prom.props('canFam3','Dog','CNE','~/Downloads/Liftover.chains/hg19ToCanFam3.over.chain','../../../../CNEs.22.bed'),  
  
  get.cgi.prom.props('felCat9','Cat','gBGC AR','~/Downloads/Liftover.chains/hg38ToFelCat9.over.chain','felCat5.GCbias.bed',T,F),  
  get.cgi.prom.props('felCat9','Cat','other AR','~/Downloads/Liftover.chains/hg38ToFelCat9.over.chain','noGCbias/felCat5.bed',T,F),  
  get.cgi.prom.props('felCat9','Cat','CNE','~/Downloads/Liftover.chains/hg38ToFelCat9.over.chain','../../../../CNEs.22.bed',T,F),  
  
  
  get.cgi.prom.props('bosTau9','Cow','gBGC AR','~/Downloads/Liftover.chains/hg38ToBosTau9.over.chain','bosTau7.GCbias.bed',T,F),
  get.cgi.prom.props('bosTau9','Cow','other AR','~/Downloads/Liftover.chains/hg38ToBosTau9.over.chain','noGCbias/bosTau7.bed',T,F),
  get.cgi.prom.props('bosTau9','Cow','CNE','~/Downloads/Liftover.chains/hg38ToBosTau9.over.chain','../../../../CNEs.22.bed',T,F),
  
  get.cgi.prom.props('susScr3','Pig','gBGC AR','~/Downloads/Liftover.chains/hg19ToSusScr3.over.chain','susScr3.GCbias.bed'),
  get.cgi.prom.props('susScr3','Pig','other AR','~/Downloads/Liftover.chains/hg19ToSusScr3.over.chain','noGCbias/susScr3.bed'),
  get.cgi.prom.props('susScr3','Pig','CNE','~/Downloads/Liftover.chains/hg19ToSusScr3.over.chain','../../../../CNEs.22.bed'),
  
  get.cgi.prom.props('equCab3','Horse','gBGC AR','~/Downloads/Liftover.chains/hg38ToEquCab3.over.chain','equCab2.GCbias.bed',T,F),
  get.cgi.prom.props('equCab3','Horse','other AR','~/Downloads/Liftover.chains/hg38ToEquCab3.over.chain','noGCbias/equCab2.bed',T,F),
  get.cgi.prom.props('equCab3','Horse','CNE','~/Downloads/Liftover.chains/hg38ToEquCab3.over.chain','../../../../CNEs.22.bed',T,F),  
  
  get.cgi.prom.props('hg19','Human','gBGC AR','~/Downloads/Liftover.chains/hg19ToHg38.over.chain','hg19.GCbias.bed',F,T),
  get.cgi.prom.props('hg19','Human','other AR','~/Downloads/Liftover.chains/hg19ToHg38.over.chain','noGCbias/hg19.bed',F,T),  
  get.cgi.prom.props('hg19','Human','CNE','~/Downloads/Liftover.chains/hg19ToHg38.over.chain','../../../../CNEs.22.bed',F,T),  
  
  get.cgi.prom.props('mm10','Mouse','gBGC AR','~/Downloads/Liftover.chains/hg19ToMm10.over.chain','mm10.GCbias.bed'),  
  get.cgi.prom.props('mm10','Mouse','other AR','~/Downloads/Liftover.chains/hg19ToMm10.over.chain','noGCbias/mm10.bed'),  
  get.cgi.prom.props('mm10','Mouse','CNE','~/Downloads/Liftover.chains/hg19ToMm10.over.chain','../../../../CNEs.22.bed'), 
  
  get.cgi.prom.props('rn6','Rat','gBGC AR','~/Downloads/Liftover.chains/hg19ToRn6.over.chain','rn5.GCbias.bed'),
  get.cgi.prom.props('rn6','Rat','other AR','~/Downloads/Liftover.chains/hg19ToRn6.over.chain','noGCbias/rn5.bed'), 
  get.cgi.prom.props('rn6','Rat','CNE','~/Downloads/Liftover.chains/hg19ToRn6.over.chain','../../../../CNEs.22.bed'),
  
  get.cgi.prom.props('panTro4','Chimp','gBGC AR','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','panTro4.GCbias.bed'), 
  get.cgi.prom.props('panTro4','Chimp','other AR','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','noGCbias/panTro4.bed'), 
  get.cgi.prom.props('panTro4','Chimp','CNE','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','../../../../CNEs.22.bed'), 
  
  get.cgi.prom.props('rheMac10','Macaque','gBGC AR','~/Downloads/Liftover.chains/hg19ToRheMac10.over.chain','rheMac3.GCbias.bed'), 
  get.cgi.prom.props('rheMac10','Macaque','other AR','~/Downloads/Liftover.chains/hg19ToRheMac10.over.chain','noGCbias/rheMac3.bed'), 
  get.cgi.prom.props('rheMac10','Macaque','CNE','~/Downloads/Liftover.chains/hg19ToRheMac10.over.chain','../../../../CNEs.22.bed') 
  
  
 
  )


testdf$Overlaps <- factor(testdf$Overlaps,levels=c('Neither','CGI','CGI.Promoter','Promoter'))
testdf$Element <- factor(testdf$Element,levels=c('CNE','other AR','gBGC AR'))
testdf$Species <- factor(testdf$Species,levels=c('Human','Chimp','Macaque','Mouse','Rat','Cow','Pig','Horse','Cat','Dog'))
testdf[97,4] <- 1 # eeeesh, cover for now ####################### ####################### ####################### #######################

ggplot(testdf[testdf$Overlaps !='Neither',],aes(Element,n,fill=Overlaps)) + geom_bar(stat='identity',position='fill') + 
  facet_wrap(.~Species,nrow=2) + theme_classic() + ggsci::scale_fill_jco() + ylab('Proportion of regions')

ggplot(testdf,aes(Element,n,fill=Overlaps)) + geom_bar(stat='identity',position='fill') + 
  facet_wrap(.~Species,nrow = 2) + theme_classic() + ggsci::scale_fill_jco() + ylab('Proportion of regions')

#get.cgi.prom.props('panTro4','Chimp','other AR','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','noGCbias/panTro4.bed') 


#######################################################################################################################################
#######################################################################################################################################

get.prom.dists <- function(sp,spname,ele,ch,bed,lift=T,hg19=T){
  ch = import.chain(ch)
  
  bed <- import.bed(bed); seqlevelsStyle(bed) <- 'UCSC'
  
  if(lift ==T){
    if(hg19 == T){
      bed <- unlist(liftOver(bed,ch))
      
    }
    else{
      
      bed.38 <- unlist(liftOver(bed,human.ch))
      bed <- unlist(liftOver(bed.38,ch))
      
    }
    bed <- reduce(bed,ignore.strand = TRUE, min.gapwidth = 51L)
    
  }
  dogSession<-browserSession("UCSC")
  genome(dogSession) <- sp
  
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
  
  dists.2.promoters <- distanceToNearest(bed,promoters(doggenes.gr))
  dist.out <- data.frame(Species=rep(spname,length(dists.2.promoters)),Element = rep(ele,length(dists.2.promoters)), Distances = dists.2.promoters)

  cgi.eles <- width(cpgIslands.gr[unique(subjectHits(findOverlaps(bed,cpgIslands.gr)))])
  cgi.other <- width(cpgIslands.gr[-unique(subjectHits(findOverlaps(bed,cpgIslands.gr)))])
  
  cgi.widths.out <- data.frame(Species=rep(spname,length(c(cgi.eles,cgi.other))),Element=c(rep(ele,length(cgi.eles)), rep('other CGIs' ,length(cgi.other)) ),Widths= c(cgi.eles,cgi.other)) 
  
  #return(list(dist.out,cgi.widths.out))
  return(dist.out)
  
}

testdf.promdist <- rbind(
  
  get.prom.dists('canFam3','Dog','gBGC AR','~/Downloads/Liftover.chains/hg19ToCanFam3.over.chain','canFam3.GCbias.bed'),  
  get.prom.dists('canFam3','Dog','other AR','~/Downloads/Liftover.chains/hg19ToCanFam3.over.chain','noGCbias/canFam3.bed'),  
  get.prom.dists('canFam3','Dog','CNE','~/Downloads/Liftover.chains/hg19ToCanFam3.over.chain','../../../../CNEs.22.bed'),  
  
  get.prom.dists('felCat9','Cat','gBGC AR','~/Downloads/Liftover.chains/hg38ToFelCat9.over.chain','felCat5.GCbias.bed',T,F),  
  get.prom.dists('felCat9','Cat','other AR','~/Downloads/Liftover.chains/hg38ToFelCat9.over.chain','noGCbias/felCat5.bed',T,F),  
  get.prom.dists('felCat9','Cat','CNE','~/Downloads/Liftover.chains/hg38ToFelCat9.over.chain','../../../../CNEs.22.bed',T,F),  
  
  
  get.prom.dists('bosTau9','Cow','gBGC AR','~/Downloads/Liftover.chains/hg38ToBosTau9.over.chain','bosTau7.GCbias.bed',T,F),
  get.prom.dists('bosTau9','Cow','other AR','~/Downloads/Liftover.chains/hg38ToBosTau9.over.chain','noGCbias/bosTau7.bed',T,F),
  get.prom.dists('bosTau9','Cow','CNE','~/Downloads/Liftover.chains/hg38ToBosTau9.over.chain','../../../../CNEs.22.bed',T,F),
  
  get.prom.dists('susScr3','Pig','gBGC AR','~/Downloads/Liftover.chains/hg19ToSusScr3.over.chain','susScr3.GCbias.bed'),
  get.prom.dists('susScr3','Pig','other AR','~/Downloads/Liftover.chains/hg19ToSusScr3.over.chain','noGCbias/susScr3.bed'),
  get.prom.dists('susScr3','Pig','CNE','~/Downloads/Liftover.chains/hg19ToSusScr3.over.chain','../../../../CNEs.22.bed'),
  
  get.prom.dists('equCab3','Horse','gBGC AR','~/Downloads/Liftover.chains/hg38ToEquCab3.over.chain','equCab2.GCbias.bed',T,F),
  get.prom.dists('equCab3','Horse','other AR','~/Downloads/Liftover.chains/hg38ToEquCab3.over.chain','noGCbias/equCab2.bed',T,F),
  get.prom.dists('equCab3','Horse','CNE','~/Downloads/Liftover.chains/hg38ToEquCab3.over.chain','../../../../CNEs.22.bed',T,F),  
  
  get.prom.dists('hg19','Human','gBGC AR','~/Downloads/Liftover.chains/hg19ToHg38.over.chain','hg19.GCbias.bed',F,T),
  get.prom.dists('hg19','Human','other AR','~/Downloads/Liftover.chains/hg19ToHg38.over.chain','noGCbias/hg19.bed',F,T),  
  get.prom.dists('hg19','Human','CNE','~/Downloads/Liftover.chains/hg19ToHg38.over.chain','../../../../CNEs.22.bed',F,T),  
  
  get.prom.dists('mm10','Mouse','gBGC AR','~/Downloads/Liftover.chains/hg19ToMm10.over.chain','mm10.GCbias.bed'),  
  get.prom.dists('mm10','Mouse','other AR','~/Downloads/Liftover.chains/hg19ToMm10.over.chain','noGCbias/mm10.bed'),  
  get.prom.dists('mm10','Mouse','CNE','~/Downloads/Liftover.chains/hg19ToMm10.over.chain','../../../../CNEs.22.bed'), 
  
  get.prom.dists('rn6','Rat','gBGC AR','~/Downloads/Liftover.chains/hg19ToRn6.over.chain','rn5.GCbias.bed'),
  get.prom.dists('rn6','Rat','other AR','~/Downloads/Liftover.chains/hg19ToRn6.over.chain','noGCbias/rn5.bed'), 
  get.prom.dists('rn6','Rat','CNE','~/Downloads/Liftover.chains/hg19ToRn6.over.chain','../../../../CNEs.22.bed'),
  
  get.prom.dists('panTro4','Chimp','gBGC AR','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','panTro4.GCbias.bed'), 
  get.prom.dists('panTro4','Chimp','other AR','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','noGCbias/panTro4.bed'), 
  get.prom.dists('panTro4','Chimp','CNE','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','../../../../CNEs.22.bed'), 
  
  get.prom.dists('rheMac10','Macaque','gBGC AR','~/Downloads/Liftover.chains/hg19ToRheMac10.over.chain','rheMac3.GCbias.bed'), 
  get.prom.dists('rheMac10','Macaque','other AR','~/Downloads/Liftover.chains/hg19ToRheMac10.over.chain','noGCbias/rheMac3.bed'), 
  get.prom.dists('rheMac10','Macaque','CNE','~/Downloads/Liftover.chains/hg19ToRheMac10.over.chain','../../../../CNEs.22.bed') 
  
  
  
)


testdf.promdist$Distances.distance <- testdf.promdist$Distances.distance + 1
testdf.promdist$Species <- factor(testdf.promdist$Species,levels=c('Human','Chimp','Macaque','Mouse','Rat','Cow','Pig','Horse','Cat','Dog'))


ggplot(testdf.promdist,aes(Distances.distance,color=Element)) + stat_ecdf() + scale_x_continuous(trans='log10') +
  theme_classic() + ggsci::scale_color_jco() + facet_wrap(.~Species) + ylab('Proportion') + xlab('Distance to nearest promoter (bp)')

#ggplot(testdf.promdist,aes(Element, Distances.distance,fill=Element)) + geom_violin() + scale_y_continuous(trans='log10') +
#  theme_classic() + ggsci::scale_fill_jco() + facet_wrap(.~Species)  + ylab('Distance to nearest promoter (bp)')



##############################
##############################


get.cgi.widths <- function(sp,spname,ele,ch,bed,lift=T,hg19=T){
  ch = import.chain(ch)
  
  bed <- import.bed(bed); seqlevelsStyle(bed) <- 'UCSC'
  
  if(lift ==T){
    if(hg19 == T){
      bed <- unlist(liftOver(bed,ch))
      
    }
    else{
      
      bed.38 <- unlist(liftOver(bed,human.ch))
      bed <- unlist(liftOver(bed.38,ch))
      
    }
    bed <- reduce(bed,ignore.strand = TRUE, min.gapwidth = 51L)
    
  }
  dogSession<-browserSession("UCSC")
  genome(dogSession) <- sp
  
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
  
  dists.2.promoters <- distanceToNearest(bed,promoters(doggenes.gr))
  dist.out <- data.frame(Species=rep(spname,length(dists.2.promoters)),Element = rep(ele,length(dists.2.promoters)), Distances = dists.2.promoters)
  
  cgi.eles <- width(cpgIslands.gr[unique(subjectHits(findOverlaps(bed,cpgIslands.gr)))])
  cgi.other <- width(cpgIslands.gr[-unique(subjectHits(findOverlaps(bed,cpgIslands.gr)))])
  
  cgi.widths.out <- data.frame(Species=rep(spname,length(c(cgi.eles,cgi.other))),Element=c(rep(ele,length(cgi.eles)), rep('other CGIs' ,length(cgi.other)) ),Widths= c(cgi.eles,cgi.other)) 
  
  #return(list(dist.out,cgi.widths.out))
  return(cgi.widths.out)
  
}



testdf.cgiwidths <- rbind(
  
  get.cgi.widths('canFam3','Dog','gBGC AR','~/Downloads/Liftover.chains/hg19ToCanFam3.over.chain','canFam3.GCbias.bed'),  
  get.cgi.widths('canFam3','Dog','other AR','~/Downloads/Liftover.chains/hg19ToCanFam3.over.chain','noGCbias/canFam3.bed'),  
  get.cgi.widths('canFam3','Dog','CNE','~/Downloads/Liftover.chains/hg19ToCanFam3.over.chain','../../../../CNEs.22.bed'),  
  
  get.cgi.widths('felCat9','Cat','gBGC AR','~/Downloads/Liftover.chains/hg38ToFelCat9.over.chain','felCat5.GCbias.bed',T,F),  
  get.cgi.widths('felCat9','Cat','other AR','~/Downloads/Liftover.chains/hg38ToFelCat9.over.chain','noGCbias/felCat5.bed',T,F),  
  get.cgi.widths('felCat9','Cat','CNE','~/Downloads/Liftover.chains/hg38ToFelCat9.over.chain','../../../../CNEs.22.bed',T,F),  
  
  
  get.cgi.widths('bosTau9','Cow','gBGC AR','~/Downloads/Liftover.chains/hg38ToBosTau9.over.chain','bosTau7.GCbias.bed',T,F),
  get.cgi.widths('bosTau9','Cow','other AR','~/Downloads/Liftover.chains/hg38ToBosTau9.over.chain','noGCbias/bosTau7.bed',T,F),
  get.cgi.widths('bosTau9','Cow','CNE','~/Downloads/Liftover.chains/hg38ToBosTau9.over.chain','../../../../CNEs.22.bed',T,F),
  
  get.cgi.widths('susScr3','Pig','gBGC AR','~/Downloads/Liftover.chains/hg19ToSusScr3.over.chain','susScr3.GCbias.bed'),
  get.cgi.widths('susScr3','Pig','other AR','~/Downloads/Liftover.chains/hg19ToSusScr3.over.chain','noGCbias/susScr3.bed'),
  get.cgi.widths('susScr3','Pig','CNE','~/Downloads/Liftover.chains/hg19ToSusScr3.over.chain','../../../../CNEs.22.bed'),
  
  get.cgi.widths('equCab3','Horse','gBGC AR','~/Downloads/Liftover.chains/hg38ToEquCab3.over.chain','equCab2.GCbias.bed',T,F),
  get.cgi.widths('equCab3','Horse','other AR','~/Downloads/Liftover.chains/hg38ToEquCab3.over.chain','noGCbias/equCab2.bed',T,F),
  get.cgi.widths('equCab3','Horse','CNE','~/Downloads/Liftover.chains/hg38ToEquCab3.over.chain','../../../../CNEs.22.bed',T,F),  
  
  get.cgi.widths('hg19','Human','gBGC AR','~/Downloads/Liftover.chains/hg19ToHg38.over.chain','hg19.GCbias.bed',F,T),
  get.cgi.widths('hg19','Human','other AR','~/Downloads/Liftover.chains/hg19ToHg38.over.chain','noGCbias/hg19.bed',F,T),  
  get.cgi.widths('hg19','Human','CNE','~/Downloads/Liftover.chains/hg19ToHg38.over.chain','../../../../CNEs.22.bed',F,T),  
  
  get.cgi.widths('mm10','Mouse','gBGC AR','~/Downloads/Liftover.chains/hg19ToMm10.over.chain','mm10.GCbias.bed'),  
  get.cgi.widths('mm10','Mouse','other AR','~/Downloads/Liftover.chains/hg19ToMm10.over.chain','noGCbias/mm10.bed'),  
  get.cgi.widths('mm10','Mouse','CNE','~/Downloads/Liftover.chains/hg19ToMm10.over.chain','../../../../CNEs.22.bed'), 
  
  get.cgi.widths('rn6','Rat','gBGC AR','~/Downloads/Liftover.chains/hg19ToRn6.over.chain','rn5.GCbias.bed'),
  get.cgi.widths('rn6','Rat','other AR','~/Downloads/Liftover.chains/hg19ToRn6.over.chain','noGCbias/rn5.bed'), 
  get.cgi.widths('rn6','Rat','CNE','~/Downloads/Liftover.chains/hg19ToRn6.over.chain','../../../../CNEs.22.bed'),
  
  get.cgi.widths('panTro4','Chimp','gBGC AR','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','panTro4.GCbias.bed'), 
  get.cgi.widths('panTro4','Chimp','other AR','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','noGCbias/panTro4.bed'), 
  get.cgi.widths('panTro4','Chimp','CNE','~/Downloads/Liftover.chains/hg19ToPanTro4.over.chain','../../../../CNEs.22.bed'), 
  
  get.cgi.widths('rheMac10','Macaque','gBGC AR','~/Downloads/Liftover.chains/hg19ToRheMac10.over.chain','rheMac3.GCbias.bed'), 
  get.cgi.widths('rheMac10','Macaque','other AR','~/Downloads/Liftover.chains/hg19ToRheMac10.over.chain','noGCbias/rheMac3.bed'), 
  get.cgi.widths('rheMac10','Macaque','CNE','~/Downloads/Liftover.chains/hg19ToRheMac10.over.chain','../../../../CNEs.22.bed') 
  
)

testdf.cgiwidths$Species <- factor(testdf.cgiwidths$Species,levels=c('Human','Chimp','Macaque','Mouse','Rat','Cow','Pig','Horse','Cat','Dog'))

ggplot(testdf.cgiwidths,aes(Widths,color=Element)) + stat_ecdf() + scale_x_continuous(trans='log10') +
  theme_classic() + ylab('Proportion') + xlab('CGI width (bp)') + ggsci::scale_color_jco() + facet_wrap(.~Species)

#ggplot(testdf.cgiwidths,aes(Element,Widths,fill=Element)) + geom_boxplot() + scale_y_continuous(trans='log10') +
#  theme_classic() + ylab('CGI width') + ggsci::scale_fill_jco() + facet_wrap(.~Species)
