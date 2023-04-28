### A script for taking cpg islands and promoters for 2 species, converting to human coordinates
## and asking whether species specific islands/promoters are more likely to overlap ars/bgc regions than cnes

##consideration- only keep genes that are refseq curated/ucsc genes in human?
library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
library(regioneR)

x <- read.table('~/Documents/Projects/Accelerated_Regions/hg19.canFam3.98.50.processed.threshold.grbs.bed',col.names = c('seqnames','start','end','loci','score','x','start2','end2','color'),skip = 1)
#x <- read.table('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/grb_boundaries_genes.hg19.txt',col.names = c('seqnames','start','end','target_gene'))
grbs.gr <- unique(makeGRangesFromDataFrame(x,keep.extra.columns = T))
grbs.gr <- grbs.gr[width(grbs.gr) > 100000]


setwd('~/Downloads/Liftover.chains.tohuman/')
human.ch <- import.chain('hg38ToHg19.over.chain')
comp.bgc.cpg <- function(ch,sp,hg38=F){
  
  Session<-browserSession("UCSC")
  genome(Session) <- sp
  
  cpgIslands<-getTable(ucscTableQuery(Session, track="cpgIslandExt", table="cpgIslandExt"))
  cpgIslands.gr<-GRanges(seqnames=cpgIslands$chrom, 
                         ranges=IRanges(start=cpgIslands$chromStart,
                                        end=cpgIslands$chromEnd),
                         gcNum=cpgIslands$gcNum,
                         cpgNum=cpgIslands$cpgNum)
  
  genes <- getTable(ucscTableQuery(Session, track="NCBI RefSeq", table="ncbiRefSeq"))
  genes.gr <- promoters(unique(GRanges(seqnames=genes$chrom, 
                                ranges=IRanges(start=genes$txStart,
                                               end=genes$txEnd),
                                name=genes$name2,
  )))
  
  
  #
  if(sp != 'hg19'){
    ch = import.chain(ch)
    
  
    genes.hg <- unlist(liftOver(genes.gr,ch))
    cpgs.hg <- unlist(liftOver(cpgIslands.gr,ch))
  
    if(hg38==T){
      genes.hg <- unlist(liftOver(genes.hg,human.ch))
      cpgs.hg <- unlist(liftOver(cpgs.hg,human.ch))
    }
        
    genes.hg <- GenomicRanges::reduce(genes.hg,ignore.strand = TRUE, min.gapwidth = 51L)
    cpgs.hg <- GenomicRanges::reduce(cpgs.hg,ignore.strand = TRUE, min.gapwidth = 51L)
  
  return(list(genes=genes.hg,cpgisland=cpgs.hg))}
  else{
    return(list(genes=genes.gr,cpgisland=cpgIslands.gr))
  }
}

#####

dog_coords <- comp.bgc.cpg('canFam3ToHg19.over.chain','canFam3',hg38 = F)
cat_coords <- comp.bgc.cpg('felCat9ToHg38.over.chain','felCat9',hg38 = T)
mouse_coords <- comp.bgc.cpg('mm10ToHg19.over.chain','mm10',hg38 = F)
rat_coords <- comp.bgc.cpg('rn6ToHg19.over.chain','rn6',hg38 = F)

human_coords <- comp.bgc.cpg('hg38ToHg19.over.chain','hg19',hg38 = F)
chimp_coords <- comp.bgc.cpg('panTro4ToHg19.over.chain','panTro4',hg38 = F)

horse_coords <- comp.bgc.cpg('equCab3ToHg38.over.chain','equCab3',hg38 = T) #doesn't lift over to hg19 blegh
cow_coords <- comp.bgc.cpg('bosTau7ToHg19.over.chain','bosTau7',hg38 = F) 

pig_coords <- comp.bgc.cpg('susScr3ToHg19.over.chain','susScr3',hg38 = F) 

cow_coords <- comp.bgc.cpg('bosTau7ToHg19.over.chain','bosTau7',hg38 = F) 


####
setwd('~/Documents/Projects/Accelerated_Regions/generalstuff/July Stuff/Accelerated_Regions.July2017/GC_bias/')
cnes.gr <- import.bed('~/Documents/Projects/Accelerated_Regions/CNEs.July17.bed')

dog.accels <- import.bed('noGCbias/canFam3.bed')
dog.bgc <- import.bed('canFam3.GCbias.bed')

cat.accels <- import.bed('noGCbias/felCat5.bed')
cat.bgc <- import.bed('felCat5.GCbias.bed')

mouse.accels <- import.bed('noGCbias/mm10.bed')
mouse.bgc <- import.bed('mm10.GCbias.bed')

rat.accels <- import.bed('noGCbias/rn5.bed')
rat.bgc  <- import.bed('rn5.GCbias.bed')

hg.accels <- import.bed('noGCbias/hg19.bed')
hg.bgc  <- import.bed('hg19.GCbias.bed')

chimp.accels <- import.bed('noGCbias/panTro4.bed')
chimp.bgc  <- import.bed('panTro4.GCbias.bed')

horse.accels <- import.bed('noGCbias/equCab2.bed')
horse.bgc  <- import.bed('equCab2.GCbias.bed')

pig.accels <- import.bed('noGCbias/susScr3.bed')
pig.bgc  <- import.bed('susScr3.GCbias.bed')

cow.accels <- import.bed('noGCbias/bosTau7.bed')
cow.bgc  <- import.bed('bosTau7.GCbias.bed')

## are promoters/cpgs more likely to overlap bgc than other accs. not quite
# is the ratio of cpgs/promoters overlapping bgc:cnes in one species dif from...

#remember to compare against cnes as well as accs
compare.bgc.cpg <- function(spA.coords,spB.coords,spA.bgc,spB.bgc,spA.acc,spB.acc){
   
   # i need to check why i did this one ^^ 
   ft1 <- fisher.test(matrix(c(numOverlaps(spA.acc,spA.coords,count.once = T),
                       numOverlaps(spA.bgc,spA.coords,count.once = T),
                       numOverlaps(spB.acc,spA.coords,count.once = T),
                       numOverlaps(spB.bgc,spA.coords,count.once = T)),ncol=2))
  
  # are species-specific promoters/cpg islands more likely to have undergone bgc than accelertion or be stable?
  #lower or means bgc is more species species
  #ie are these 'new' promoters/cpg islands 
  shared <- spA.coords[unique(queryHits(findOverlaps(spA.coords,spB.coords)))]
  spA_only <- spA.coords[-unique(queryHits(findOverlaps(spA.coords,spB.coords)))]
  
  ft2 <- fisher.test(matrix(c(
    numOverlaps(shared,spA.bgc),
    numOverlaps(spA_only,spA.bgc),
    numOverlaps(shared,spA.acc),
    numOverlaps(spA_only,spA.acc)
  ),ncol=2))
  #
  ###### 
  # are species-specific regions (cpg islands or promoters) more likely to be in grbs than conserved ones...
  #
  ft3 <- fisher.test(matrix(c(
    length(spA_only[unique(queryHits(findOverlaps(spA_only,grbs.gr)))]),
    length(spA_only[-unique(queryHits(findOverlaps(spA_only,grbs.gr)))]),
    length(shared[unique(queryHits(findOverlaps(shared,grbs.gr)))]),
    length(shared[-unique(queryHits(findOverlaps(shared,grbs.gr)))])
  ),ncol=2))
  
  #### are the ratios of species-specific promoters/cpg islands overlapping bgc ars different to species-shared overlapping...
  #how is this different from 2? This one doesn't compare to CNEs- it's comparing to other eg cpgislands
  spA_only.bcg <- length(spA_only[unique(queryHits(findOverlaps(spA_only,spA.bgc)))])
  spA_only.other <- length(spA_only[-unique(queryHits(findOverlaps(spA_only,spA.bgc)))])
  
  shared.bcg <- length(shared[unique(queryHits(findOverlaps(shared,spA.bgc)))])
  shared.other <- length(shared[-unique(queryHits(findOverlaps(shared,spA.bgc)))])
  
  ft4 <- fisher.test(matrix(c(spA_only.bcg,spA_only.other,
                       shared.bcg,shared.other),ncol=2))
  return(list(ft1,ft2,ft3,ft4))
}

# output 2 and 4 most interesting i think
compare.bgc.cpg(dog_coords[[1]],cat_coords[[1]],dog.bgc,cat.bgc,cnes.gr,cnes.gr) #dog.accels,cat.accels)
compare.bgc.cpg(cat_coords[[1]],dog_coords[[1]],cat.bgc,dog.bgc,cnes.gr,cnes.gr)

compare.bgc.cpg(dog_coords[[2]],cat_coords[[2]],dog.bgc,cat.bgc,cnes.gr,cnes.gr)
compare.bgc.cpg(cat_coords[[2]],dog_coords[[2]],cat.bgc,dog.bgc,cnes.gr,cnes.gr)
#
#sooo... both cats and dogs are gaining a lot of cpg islands?
#but for promoters.. no difference with cat, whilst dog has way more promoters linked that are species specific
#but aren't caused by gc bias?

compare.bgc.cpg(horse_coords[[2]],dog_coords[[2]],horse.bgc,dog.bgc,cnes.gr,cnes.gr)

suppressWarnings(compare.bgc.cpg2(cow_coords[[2]],pig_coords[[2]],cow.bgc,pig.bgc,cnes.gr,cnes.gr)
)
##
compare.bgc.cpg(mouse_coords[[1]],rat_coords[[1]],mouse.bgc,rat.bgc,cnes.gr,cnes.gr)
compare.bgc.cpg(rat_coords[[1]],mouse_coords[[1]],rat.bgc,mouse.bgc,cnes.gr,cnes.gr)

compare.bgc.cpg(mouse_coords[[2]],rat_coords[[2]],mouse.bgc,rat.bgc,cnes.gr,cnes.gr)
compare.bgc.cpg(rat_coords[[2]],mouse_coords[[2]],rat.bgc,mouse.bgc,cnes.gr,cnes.gr)

#

compare.bgc.cpg(human_coords[[1]],chimp_coords[[1]],hg.bgc,chimp.bgc,cnes.gr,cnes.gr)
compare.bgc.cpg(chimp_coords[[1]],human_coords[[1]],chimp.bgc,hg.bgc,cnes.gr,cnes.gr)

compare.bgc.cpg(human_coords[[2]],chimp_coords[[2]],hg.bgc,chimp.bgc,cnes.gr,cnes.gr)
compare.bgc.cpg(chimp_coords[[2]],human_coords[[2]],chimp.bgc,hg.bgc,cnes.gr,cnes.gr)

suppressWarnings(compare.bgc.cpg(pig_coords[[2]],cat_coords[[2]],pig.bgc,cat.bgc,cnes.gr,cnes.gr))

#use cnes rather than accels
#most species aren't gaining a bunch of cpg islands/promoters... but dogs/cats are
###################
###################
"""
fisher.test(matrix(c(numOverlaps(cat.accels,cat_coords[[1]]),
                     numOverlaps(cat.bgc,cat_coords[[1]]),
                     numOverlaps(dog.accels,cat_coords[[1]]),
                     numOverlaps(dog.bgc,cat_coords[[1]])),ncol=2))

# are species-specific promoters/cpg islands more likely to have undergone bgc than accelertion or be stable?
shared <- cat_coords[[1]][unique(queryHits(findOverlaps(cat_coords[[1]],dog_coords[[1]])))]
cat_only <- cat_coords[[1]][-unique(queryHits(findOverlaps(cat_coords[[1]],dog_coords[[1]])))]

fisher.test(matrix(c(
  numOverlaps(shared,cat.bgc),
  numOverlaps(cat_only,cat.bgc),
  numOverlaps(shared,cat.accels),
  numOverlaps(cat_only,cat.accels)
),ncol=2))
#
###### so (for cats) a higher proportion of CpG islands (~2:1) are species specific 
###### whilst the opposite is seen for CNEs (~1:2)
###### it looks like similar may be the case for promoters (124:69) vs (54:43) but not big enough difference to be sure
# 

fisher.test(matrix(c(
  length(cat_only[unique(queryHits(findOverlaps(cat_only,grbs.gr)))]),
  length(cat_only[-unique(queryHits(findOverlaps(cat_only,grbs.gr)))]),
  length(shared[unique(queryHits(findOverlaps(shared,grbs.gr)))]),
  length(shared[-unique(queryHits(findOverlaps(shared,grbs.gr)))])
),ncol=2))

#for cat, species unique promoters are more likely to be in grbs. but lineage specific cpg islands are not
#stronger for genes than cpgs. both cases, new promoters/cpg islands more likely to be in grbs
#but is it stronger for bgc than not?

cat_only.bcg <- length(cat_only[unique(queryHits(findOverlaps(cat_only,cat.bgc)))])
cat_only.other <- length(cat_only[-unique(queryHits(findOverlaps(cat_only,cat.bgc)))])

shared.bcg <- length(shared[unique(queryHits(findOverlaps(shared,cat.bgc)))])
shared.other <- length(shared[-unique(queryHits(findOverlaps(shared,cat.bgc)))])

fisher.test(matrix(c(cat_only.bcg,cat_only.other,
                     shared.bcg,shared.other),ncol=2))

#nop, opposite. cat bgc promoters less likely to be found in grbs. but cpg islands more likely? 
#if use accels rather than bgc, probs switch around i think, which is nice
"""