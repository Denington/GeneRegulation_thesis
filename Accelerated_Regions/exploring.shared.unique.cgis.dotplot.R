### A script for taking cpg islands and promoters for 2 species, converting to human coordinates
## and asking whether species specific islands/promoters are more likely to overlap ars/bgc regions than cnes

##consideration- only keep genes that are refseq curated/ucsc genes in human?
library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
library(regioneR)

setwd('~/Downloads/Liftover.chains.tohuman/')

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
  
  #genes <- getTable(ucscTableQuery(Session, track="NCBI RefSeq", table="ncbiRefSeq"))
  #genes.gr <- promoters(unique(GRanges(seqnames=genes$chrom, 
  #                                     ranges=IRanges(start=genes$txStart,
  #                                                    end=genes$txEnd),
  #                                     name=genes$name2,
  #)))
  
  
  #
  if(sp != 'hg19'){
    ch = import.chain(ch)
    
    
    #genes.hg <- unlist(liftOver(genes.gr,ch))
    cpgs.hg <- unlist(liftOver(cpgIslands.gr,ch))
    
    if(hg38==T){
      #genes.hg <- unlist(liftOver(genes.hg,human.ch))
      cpgs.hg <- unlist(liftOver(cpgs.hg,human.ch))
    }
    
    #genes.hg <- GenomicRanges::reduce(genes.hg,ignore.strand = TRUE, min.gapwidth = 51L)
    cpgs.hg <- GenomicRanges::reduce(cpgs.hg,ignore.strand = TRUE, min.gapwidth = 51L)
    
    return(list(genes='genes.hg',cpgisland=cpgs.hg))}
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
macaque_coords <- comp.bgc.cpg('rheMac3ToHg19.over.chain','rheMac3',hg38 = F)

horse_coords <- comp.bgc.cpg('equCab3ToHg38.over.chain','equCab3',hg38 = T) #doesn't lift over to hg19 blegh
cow_coords <- comp.bgc.cpg('bosTau7ToHg19.over.chain','bosTau7',hg38 = F) 

pig_coords <- comp.bgc.cpg('susScr3ToHg19.over.chain','susScr3',hg38 = F) 

#dolphin_coords <- comp.bgc.cpg('turTru2ToHg38.over.chain','susScr3',hg38 = T) 
#bat_coords <- comp.bgc.cpg('myoLuc2ToHg38.over.chain','susScr3',hg38 = T) 

#g.pig_coords <- comp.bgc.cpg('cavPor3ToHg38.over.chain','cavPor3',hg38 = T) 

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

macaque.accels <- import.bed('noGCbias/rheMac3.bed')
macaque.bgc  <- import.bed('rheMac3.GCbias.bed')

horse.accels <- import.bed('noGCbias/equCab2.bed')
horse.bgc  <- import.bed('equCab2.GCbias.bed')

pig.accels <- import.bed('noGCbias/susScr3.bed')
pig.bgc  <- import.bed('susScr3.GCbias.bed')

cow.accels <- import.bed('noGCbias/bosTau7.bed')
cow.bgc  <- import.bed('bosTau7.GCbias.bed')

## are promoters/cpgs more likely to overlap bgc than other accs. not quite
# is the ratio of cpgs/promoters overlapping bgc:cnes in one species dif from...

#remember to compare against cnes as well as accs
compare.bgc.cpg <- function(species,spA.coords,spB.coords,spA.bgc,spA.acc,spA.CNEs){
  

  # are species-specific promoters/cpg islands more likely to have undergone bgc than accelertion or be stable?
  #lower or means bgc is more species species
  #ie are these 'new' promoters/cpg islands 
  shared <- spA.coords[unique(queryHits(findOverlaps(spA.coords,spB.coords)))]
  spA_only <- spA.coords[-unique(queryHits(findOverlaps(spA.coords,spB.coords)))]
  
  bgc.ls.spratio <- 100*numOverlaps(spA_only,spA.bgc)/(numOverlaps(spA_only,spA.bgc)+ numOverlaps(shared,spA.bgc))
  acc.ls.spratio <- 100*numOverlaps(spA_only,spA.acc)/(numOverlaps(spA_only,spA.acc)+ numOverlaps(shared,spA.acc))
  cne.ls.spratio <- 100*numOverlaps(spA_only,spA.CNEs)/(numOverlaps(spA_only,spA.CNEs)+ numOverlaps(shared,spA.CNEs))
  
  spA_only.other <- length(spA_only[-unique(queryHits(findOverlaps(spA_only,spA.CNEs)))])
  shared.other <- length(shared[-unique(queryHits(findOverlaps(shared,spA.CNEs)))])
  other.ratio <- 100*spA_only.other /(spA_only.other + shared.other)

  return(data.frame(Species=species,Element=c('gBGC AR','other AR','CNE','not CNEs'),Ratios=c(bgc.ls.spratio,acc.ls.spratio,cne.ls.spratio,other.ratio)))
}

# output 2 and 4 most interesting i think
bgc.r.df <- rbind(
  suppressWarnings(compare.bgc.cpg('Dog',dog_coords[[2]],cat_coords[[2]],dog.bgc,dog.accels,cnes.gr)),
  suppressWarnings(compare.bgc.cpg('Cat',cat_coords[[2]],horse_coords[[2]],cat.bgc,cat.accels,cnes.gr)),
  suppressWarnings(compare.bgc.cpg('Horse',horse_coords[[2]],dog_coords[[2]],horse.bgc,horse.accels,cnes.gr)),
  suppressWarnings(compare.bgc.cpg('Pig',pig_coords[[2]],horse_coords[[2]],pig.bgc,pig.accels,cnes.gr)),
  suppressWarnings(compare.bgc.cpg('Cow',cow_coords[[2]],dog_coords[[2]],cow.bgc,cow.accels,cnes.gr)),
  suppressWarnings(compare.bgc.cpg('Cow',cow_coords[[2]],dog_coords[[2]],cow.bgc,cow.accels,cnes.gr)),
  
)

#suppressWarnings(compare.bgc.cpg('Chimp',chimp_coords[[2]],human_coords[[2]],chimp.bgc,chimp.accels,cnes.gr))
#suppressWarnings(compare.bgc.cpg('Human',human_coords[[2]],chimp_coords[[2]],hg.bgc,hg.accels,cnes.gr))
#suppressWarnings(compare.bgc.cpg('Macaque',macaque_coords[[2]],human_coords[[2]],macaque.bgc,macaque.accels,cnes.gr))
#suppressWarnings(compare.bgc.cpg('Mouse',mouse_coords[[2]],rat_coords[[2]],mouse.bgc,mouse.accels,cnes.gr))
#suppressWarnings(compare.bgc.cpg('Rat',rat_coords[[2]],mouse_coords[[2]],rat.bgc,rat.accels,cnes.gr))

bgc.r.df$Element <- factor(bgc.r.df$Element,levels=c('not conserved','gBGC AR','other AR','CNE'))
bgc.r.df$Species <- factor(bgc.r.df$Species,levels=c('Cow','Pig','Horse','Cat','Dog'))

ggplot(bgc.r.df,aes(Species,Ratios,color=Element)) + geom_point(size=4) + theme_bw() + ggsci::scale_color_jco() +
  ylab('% of CGIs that are Lineage Specific') + theme(text = element_text(size = 15)) 
ggplot(bgc.r.df,aes(Species,Ratios,fill=Element)) + geom_bar(stat='identity',position='dodge') + theme_bw() + ggsci::scale_fill_jco() +
  ylab('% of CGIs that are Lineage Specific') + theme(text = element_text(size = 15)) 

#suppressWarnings(compare.bgc.cpg('Dog',dog_coords[[2]],cat_coords[[2]],dog.bgc,dog.accels,cnes.gr)),
#suppressWarnings(compare.bgc.cpg('Cat',cat_coords[[2]],dog_coords[[2]],cat.bgc,cat.accels,cnes.gr)),
#suppressWarnings(compare.bgc.cpg('Horse',horse_coords[[2]],dog_coords[[2]],horse.bgc,horse.accels,cnes.gr)),
#suppressWarnings(compare.bgc.cpg('Human',human_coords[[2]],chimp_coords[[2]],hg.bgc,hg.accels,cnes.gr)),
#suppressWarnings(compare.bgc.cpg('Chimp',chimp_coords[[2]],human_coords[[2]],chimp.bgc,chimp.accels,cnes.gr)),
#suppressWarnings(compare.bgc.cpg('Pig',pig_coords[[2]],horse_coords[[2]],pig.bgc,pig.accels,cnes.gr))

numOverlaps(rat.bgc,rat_coords[[2]])