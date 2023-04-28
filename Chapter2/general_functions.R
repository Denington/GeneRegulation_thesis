#General script for making heatmaps
#also getting beds/ bedpes/ granges in an appropriate for doing so
#Just wrapping old scripts here
#loads conds.gr.list when called
#ok so... 

library(genomation)
library(GenomicRanges)
library(rtracklayer)



######################################################################################
######################################################################################
# Getting gene expression for conditions, linking to locations

# saves in a list called 'conds.gr.list'
# condtions saved as 'dev' 'dko' 'ctcf_ko' and 'rad21_ko'

#originally (and unchanged from) 'ctcf_rad21_at_genes.r'
######################################################################################
######################################################################################


load('~/Documents/Projects/Thymocyte_HiC/Thesis/scripts/mm9_gene_expression.rdata.RData')
#take wtwt,wtko etc from insulation_differences_gene_expression

difexp_to_granges <- function(genechanges){
  changes <- genechanges #[wtwt$padj < 0.05,]
  changes.gr <- mm9_pc_genes.gr[mm9_pc_genes.gr$ens_id %in% genechanges$ensembl_gene_id ]
  changes.gr <- changes.gr[!duplicated(changes.gr$ens_id)]
  changes.gr <- promoters(changes.gr,upstream = 100,downstream = 100 ) #YES!
  
  changes <- changes[changes$ensembl_gene_id %in% changes.gr$ens_id,]
  changes <- changes[!duplicated(changes$ensembl_gene_id),]
  
  changes <- changes[order(changes$ensembl_gene_id),]
  changes.gr <- changes.gr[order(changes.gr$ens_id)]
  changes.gr <- promoters(changes.gr,upstream = 100,downstream = 100)
  
  changes.gr$log2FoldChange <- changes$log2FoldChange
  changes.gr$padj  <- changes$padj
  changes.gr$baseMean  <- changes$baseMean
  
  return(changes.gr)
}
#should also get a 'disreg in both' to see if the marks are particularly disreged at genes disreged in >1 cond
setwd('~/Documents/Projects/Thymocyte_HiC/DEseq2_in_vivo_Development_YaRNAseq/')

wtko <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_DKO_CD69nDPWT1-3_vs_CD69nDPDKO1-3.csv',header=T)
wtko <- wtko[!is.na(wtko$padj),]

wtwt <- read.csv('04_CD69nDPWT1-3_vs_CD69pCD4SPWT1-2.csv',header=T)
wtwt <- wtwt[!is.na(wtwt$padj),] 

setwd('~/Documents/Projects/Thymocyte_HiC/DEseq2_in_vivo_Development_YaRNAseq/')

wtctcf <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_CTCFKO_CD69nDPWT4-5_vs_CD69nDPCTCFKO1-2.csv',header=T)
wtctcf <- wtctcf[!is.na(wtctcf$padj),]

wtrad21 <- read.csv('../DEseq2_YaRNAseq/DEseq2_Output_Rad21KO_DPWT1-2_vs_DPRad21KO1-2.csv',header=T)
wtrad21 <- wtrad21[!is.na(wtrad21$padj),] 

conds <- list(dev=wtwt,dko=wtko,ctcf_ko=wtctcf,rad21_ko=wtrad21)
conds.gr.list <- lapply(conds,difexp_to_granges)

###

#################################################################################
#################################################################################
# functions for setting up bedpe files for heatmaps
# then either centering on edge or not
# returns granges
#################################################################################
#################################################################################


read_bedpe <- function(bedpe,header){
  if(header == T){
    x <- read.table(bedpe,header = T)
    x$chrom <- paste0('chr',x$chr1)
    y <- makeGRangesFromDataFrame(x,seqnames.field = 'chrom',start.field = 'x1',end.field = 'x2',keep.extra.columns = T)}
  else{
    x <- read.table(bedpe,header = F)
    x$chrom <- paste0('chr',x$V1)
    y <- makeGRangesFromDataFrame(x,seqnames.field = 'chrom',start.field = 'V2',end.field = 'V3',keep.extra.columns = T)
    
  }
  return(y)
}


read_bedpe.split.start_end <- function(bedpe,header){
  if(header == T){
    x <- read.table(bedpe,header = T)
    x$chrom <- paste0('chr',x$chr1)
    x.starts <- makeGRangesFromDataFrame(x,seqnames.field = 'chrom',start.field = 'x1',end.field = 'x2')
    end(x.starts) <- start(x.starts) + 10 #should it be one bin? hmm or centered 
    
    x.ends <- makeGRangesFromDataFrame(x,seqnames.field = 'chrom',start.field = 'x1',end.field = 'x2')
    start(x.ends) <- end(x.ends) 
    end(x.ends) <- start(x.ends) + 10}
  else{
    x <- read.table(bedpe,header = F)
    x$chrom <- paste0('chr',x$V1)
    x.starts <- makeGRangesFromDataFrame(x,seqnames.field = 'chrom',start.field = 'V2',end.field = 'V3')
    end(x.starts) <- start(x.starts) + 10 #should it be one bin? hmm or centered 
    
    x.ends <- makeGRangesFromDataFrame(x,seqnames.field = 'chrom',start.field = 'V2',end.field = 'V3')
    start(x.ends) <- end(x.ends) 
    end(x.ends) <- start(x.ends) + 10
  }
  return(c(x.starts,x.ends))
}

######################################################################################
######################################################################################
# take granges of e.g. subcompartments/ tads, w/e
# makes funnel or boundary based
# get scoreMatrix for heatmatrixing 
######################################################################################
######################################################################################


get_ScoreMatrix.Funnel <- function(domains,mark_file,size,windows,mark_format){
  domains.resize <- resize(domains[order(width(domains), decreasing=T)] ,size,'center')
  sm <- ScoreMatrixBin(mark_file,domains.resize,bin.num = windows)
  if( mark_format == 'bw' | mark_format == 'bigwig' | mark_format == 'BigWig'){
    sm <- ScoreMatrixBin(mark_file,domains.resize,bin.num = windows) }
  else{ #assume is a granges with score as weight column
    sm <- ScoreMatrixBin(mark_file,domains.resize,bin.num = windows,weight.col = 'score') }
  
  return(sm)
}

get_ScoreMatrix.boundaryCentered <- function(domains,markfile,size,win_num,mark_format){
  #don't need this if have read in bedpe with split.start.end
  domains.starts <- domains
  end(domains.starts) <- start(domains.starts) + 1 #should it be one bin? hmm or centered 
  strand(domains.starts) <- '+'
  
  domains.ends <- domains
  start(domains.ends) <- end(domains.ends) 
  end(domains.ends) <- start(domains.ends) + 1
  strand(domains.ends) <- '-'
  
  domains <- c(domains.starts,domains.ends)
  
  domain_windows <- resize(domains,width = size,fix = 'center')
  if( mark_format == 'bw' | mark_format == 'bigwig' | mark_format == 'BigWig'){
    scores.sm <- ScoreMatrixBin(markfile,domain_windows,bin.num = win_num,strand.aware = T) }
  else{ #assume is a granges with score as weight column
    scores.sm <- ScoreMatrixBin(markfile,domain_windows,bin.num = win_num,strand.aware = T,weight.col = 'score') }
  
  return(scores.sm)
}

######################################################################################
######################################################################################
# functions for subsetting gene granges 
# (conds.gr.list)
# main idea e.g. splitting contact domains into containing up or down regulated genes
######################################################################################
######################################################################################

split.by.expression <- function(genes,cond, pval=0.05,foldchange = 0 ){
  
  down.genes <-  genes[[cond]][conds.gr.list[[cond]]$padj < pval & genes[[cond]]$log2FoldChange < foldchange,]
  up.genes <-  genes[[cond]][conds.gr.list[[cond]]$padj < pval & genes[[cond]]$log2FoldChange > -foldchange,]
  stable.genes <- genes[[cond]][genes[[cond]]$padj > pval,]
  
  return(list(down.genes,up.genes,stable.genes))
}

overlap.domains.features <- function(domains, features, feature.islist=T, unique.domains=T){
  # Assumes if list of features, is in form downreg, upreg, stable. Unique domains determines whether 
  # domain can only include up or down genes
  #gives e.g. domains that contain upreg genes, downreg genes and stable genes
  # for the features.islist == T, expects split.by.expression output as feature input
  if( feature.islist == F){
  out = unique(domains[queryHits(findOverlaps(domains,features))])}
  else {
    down = unique(domains[queryHits(findOverlaps(domains,features[[1]]))])
    up = unique(domains[queryHits(findOverlaps(domains,features[[2]]))])
    stable = unique(domains[queryHits(findOverlaps(domains,features[[3]]))])
    out = list(down,up,stable)
    if(unique.domains == T){
      down.only <- down[! down %in% up ]
      up.only <- up[! up %in% down]
      
      stable.only <- stable[! (stable %in% down | stable %in% up)]
      no.expressed.genes <- domains[!(domains %in% down | domains %in% up | domains %in% stable)]

      out = list(down.only,up.only,stable.only)
      }
    }
  return(out)
}

####

overlap.domains.2features <- function(domains, features1,features2, unique.domains=T){
  # same as overlap.domains.features, but a) only for list of features, b) takes 2 different feature lists
  # Unique domains determines whether domain can only include up or down genes
  #gives e.g. domains that contain upreg genes in both sets, downreg genes in both sets etc and stable genes
  #expects split.by.expression output as feature inputs
  down1 = unique(domains[queryHits(findOverlaps(domains,features1[[1]]))])
  up1 = unique(domains[queryHits(findOverlaps(domains,features1[[2]]))])
  stable1 = unique(domains[queryHits(findOverlaps(domains,features1[[3]]))])

  down2 = unique(domains[queryHits(findOverlaps(domains,features2[[1]]))])
  up2 = unique(domains[queryHits(findOverlaps(domains,features2[[2]]))])
  stable2 = unique(domains[queryHits(findOverlaps(domains,features2[[3]]))])

  if(unique.domains == T){
    down1.only <- down1[! down1 %in% up1 ]
    up1.only <- up1[! up1 %in% down1]
    stable1.only <- stable[! (stable %in% down1 | stable %in% up1)]
    #no.expressed.genes1 <- domains[!(domains %in% down1 | domains %in% up1 | domains %in% stable1)]
    
    down2.only <- down2[! down2 %in% up2 ]
    up2.only <- up2[! up2 %in% down2]
    stable2.only <- stable[! (stable %in% down2 | stable %in% up2)]
    #no.expressed.genes2 <- domains[!(domains %in% down2 | domains %in% up2 | domains %in% stable2)]
    
    down1 <- down1.only
    up1 <- up1.only
    stable1 <- stable1.only
    down2 <- down2.only
    up2 <- up2.only
    stable2 <- stable2.only
    
    }
  up.up <- up1[queryHits(findOverlaps(up1,up2,type = 'equal'))] #i think this is how equal works?
  down.down <- down1[queryHits(findOverlaps(down1,down2,type = 'equal'))]
  
  up.down <- up1[queryHits(findOverlaps(up1,down2,type = 'equal'))]
  down.up <- down1[queryHits(findOverlaps(down1,up2,type = 'equal'))]
  out <- list(upup= up.up, downdown = down.down, updown= up.down, downup = down.up)
  return(out)
}

##


overlap.outside.domains.features <- function(domains, features, distfrom, feature.islist=T, unique.domains=T){
  # Assumes if list of features, is in form downreg, upreg, stable. Unique domains determines whether 
  # domain can only include up or down genes
  #gives e.g. domains that have upreg genes, downreg genes and stable genes <distance from the domain (not including inside)
  # for the features.islist == T, expects split.by.expression output as feature input
  if( feature.islist == F){
    inside = unique(domains[queryHits(findOverlaps(domains,features))])
    allnear = unique(domains[queryHits(findOverlaps(domains,features,maxgap = distfrom ))])
    out = allnear[-queryHits(findOverlaps(allnear,inside))] 
  }
  else {
    down.inside = unique(domains[queryHits(findOverlaps(domains,features[[1]]))])
    down.allnear = unique(domains[queryHits(findOverlaps(domains,features[[1]],maxgap = distfrom))])
    down = down.allnear[-queryHits(findOverlaps(down.allnear,down.inside))]
    
    up.inside = unique(domains[queryHits(findOverlaps(domains,features[[2]]))])
    up.allnear = unique(domains[queryHits(findOverlaps(domains,features[[2]],maxgap = distfrom))])
    up = up.allnear[-queryHits(findOverlaps(up.allnear,up.inside))]
    
    stab.inside = unique(domains[queryHits(findOverlaps(domains,features[[3]]))])
    stab.allnear = unique(domains[queryHits(findOverlaps(domains,features[[3]],maxgap = distfrom))])
    stable = stab.allnear[-queryHits(findOverlaps(stab.allnear,stab.inside))]
    
    out = list(down,up,stable)
    if(unique.domains == T){
      down.only <- down[! down %in% up ]
      up.only <- up[! up %in% down]
      
      stable.only <- stable[! (stable %in% down | stable %in% up)]
      no.expressed.genes <- domains[!(domains %in% down | domains %in% up | domains %in% stable)]
      
      out = list(down.only,up.only,stable.only)
    }
  }
  return(out)
}


####################################################################################################################
####################################################################################################################
## Adapting regioneR resampleRegions and numOverlaps to work with loops
## so can get empirical Pvals for prom-prom loop interactions etc

num.link.Olaps.loop.fts  <- function(A,B,C=list(),count.once=T, ...){
  #how many loops have a feature e.g. CTCF/ promoters on each anchor
  #requires feature at each end of loop (if B and C, B one end, C other, otherwise B both ends e.g. promoter-promoter loops)
  #loop rather than feature centric 
  if(length(C) > 0){
    o <- A[unique(linkOverlaps(A,B,C))$query]
    
  }
  else{
    o <- A[unique(linkOverlaps(A,B))$query]
  }
  if(count.once == T){
    return(length(unique(o)))}
  else{
    return(length(o))
  }
  }

num.link.Olaps.ft.loop  <- function(A,B,C = list(),count.once=T, ...){
  #how many times does a feature overlap a loop, either at both ends, or with a 2nd feature at the other end
  #requires feature at each end of loop (if B and C, B one end, C other, otherwise B both ends e.g. promoter-promoter loops)
  #feature centric (e.g. dereg genes more likely than stable to be in P-P loops?)
  #B is loop, A is feature of interest
  if(length(C) > 0){
    o <- A[linkOverlaps(B,A,C)$subject1]
    
  }
  else{
    o <- A[linkOverlaps(B,A)$subject1]
  }
  if(count.once == T){
    return(length(unique(o)))}
  else{
    return(length(o))
  }
}


resampleRegions.loops <- function(A, universe, per.chromosome=FALSE, ...) { 
  #still unhappy when per.chromosome, just doing whole genome for now..
  if(!hasArg(A)) stop("A is missing")
  if(!hasArg(universe)) stop("universe is missing")
  if(!is.logical(per.chromosome)) stop("per.chromosome must be logical")
  
  if(per.chromosome){
    chrResample <- function(chr) {
      Achr <- A[seqnames(anchors(cd69nDPCTCFKO.loops)$first) == chr & seqnames(anchors(cd69nDPCTCFKO.loops)$second) == chr]
      universe.chr <- universe[seqnames(anchors(universe)$first) == chr & seqnames(anchors(universe)$second) == chr]
      resample.chr <- universe.chr[sample(1:length(universe.chr), length(Achr))]
      return(resample.chr)
    }
    
    chr.resampled <- lapply(as.list(seqlevels(anchors(A)$first)), chrResample)
    resampled <- do.call(c, chr.resampled)
    
  }else{
    resampled <- universe[sample(1:length(universe), length(A))]
  }
  
  return(resampled)
  
}


#

permTest.loops <- function(A, ntimes=100, randomize.function, evaluate.function, alternative="auto", min.parallel=1000, force.parallel=NULL, randomize.function.name=NULL, evaluate.function.name=NULL, verbose=FALSE, ...) {
  
  #check arguments
  alternative<-match.arg(alternative,c("less","greater", "auto"))
  if(!hasArg(A)) stop("A is missing")
  if(!is.numeric(ntimes)) stop("ntime must be numeric")
  if(!hasArg(randomize.function)) stop("randomize.function is missing")
  if(!is.function(randomize.function)) stop("randomize.function must be a function")
  if(!hasArg(evaluate.function)) stop("evaluate.function is missing")
  if(!(is.function(evaluate.function) | is.list(evaluate.function))) stop("evaluate.function must be a function")
  if(!is.numeric(min.parallel)) stop("min.parallel must be numeric")
  if(ntimes<100) print(paste0("Note: The minimum p-value with only ",ntimes," permutations is ",1/(ntimes+1),". You should consider increasing the number of permutations."))
  
  #A <- toGRanges(A) #does nothing if already a GRanges object
  
  
  if(!is.null(force.parallel)) {
    doParallel <- force.parallel
  } else {
    doParallel <- (length(A)*ntimes > min.parallel)
    if(verbose) message("Auto-setting parallel computing to: ", doParallel)
  }
  
  
  #Evaluation Function: get the function name and convert to list if its not yet
  if(!is.list(evaluate.function)) { #if it's a single function
    if(is.null(evaluate.function.name)) {
      evaluate.function.name <- as.character(match.call()["evaluate.function"])
    }
    ef <- list()
    ef[[evaluate.function.name]] <-  evaluate.function
    evaluate.function <- ef
  } else { #if it's a list of functions
    if(!is.null(evaluate.function.name)) { #if names were explicitely provided
      names(evaluate.function) <- evaluate.function.name
    } else { #try to leave the current names or create new ones if no names present
      if(is.null(names(evaluate.function))) {
        names(evaluate.function) <- paste0("Function", seq_len(length(evaluate.function)))
      }
    }
  }
  
  #Randomization Function: Get a name
  if(is.null(randomize.function.name)) {
    randomize.function.name <- match.call()["randomize.function"]
  }
  
  #Start the permutation test
  #compute the evaluation function(s) using the original region set A
  original.evaluate <- sapply(seq_len(length(evaluate.function)), function(i,...) {return(evaluate.function[[i]](A,...))}, ...)
  
  if(!is.numeric(original.evaluate)) {
    stop(paste0("The evaluation function must return a numeric value but it returned an object of class ", class(original.evaluate)))
  }
  
  if(verbose) {
    #WARNING: to give some visual information about the computation done, we use a Progress Bar. However, since the GUI will not be updated
    #whilst in multi core processing (mclapply), we partition the work in chunks and repeatedly call the parallel computation, once per chunk.
    #Between chunks, the progess bar will be updated. The problem is that the total computation time might increase due to waiting and joining
    #multiple times.
    #Create the progress bar
    
    pb <- txtProgressBar(min = 0, max = ntimes, style = 3)
    setTxtProgressBar(pb, 0)
  }
  
  #define the function to create and evaluate the random sets
  randomize_and_evaluate <- function(foo, ...) {
    #randomize
    randomA <- randomize.function(A,...)
    #evaluate the random region set    
    if(verbose) {
      setTxtProgressBar(pb, foo)
    }
    
    #compute the evaluation function(s) using the RANDOMIZED region set randomA
    rand.evaluate <- sapply(seq_len(length(evaluate.function)), function(i, ...) {return(evaluate.function[[i]](randomA,...))}, ...)
    
    return(rand.evaluate)
  }
  
  #create the random sets and evaluate them  
  if(doParallel) {
    if(verbose) { #if verbose, we will do the computations in chunks and update the progress bar in between
      random.evaluate <- numeric()
      chunk.size <- max(round(ntimes/100+1), 10)
      e <- 0
      done <- FALSE
      while(!done) {
        s <- e + 1
        e <- s + chunk.size
        if(e >= ntimes) {
          e <- ntimes
          done <- TRUE
        }
        random.evaluate <- rbind(random.evaluate, do.call(rbind, mclapply(c(s:e), randomize_and_evaluate, ...)))
        setTxtProgressBar(pb, e)
      }    
    } else { #if not verbose, just do it
      random.evaluate <- do.call(rbind, mclapply(seq_len(ntimes), randomize_and_evaluate, ...))
    }
  } else {
    random.evaluate <- do.call(rbind, lapply(seq_len(ntimes), randomize_and_evaluate, ...))
  }
  
  
  #The simulation process has finished. Now build a permTestResults object for each evaluate.function
  results <- list()
  for(i in seq_len(length(evaluate.function))) {
    #Get the data for the i-th function
    func.name <- names(evaluate.function)[i]
    orig.ev <- original.evaluate[i]
    rand.ev <- random.evaluate[,i]
    
    
    #warn if any NA or NaN
    num.nas <- length(which(is.na(rand.ev)))
    num.valid.values <- ntimes-num.nas
    
    if(num.valid.values < ntimes) {
      if(num.valid.values>0) {
        warning(paste0(num.nas, " iterations returned NA or NaN. Only "  ," iterations have been used to compute the p-value."))
      } else {
        warning(paste0("All ", num.nas, " iterations returned NA or NaN. No valid values returned. It is not possible to compute the p-value nor z-score."))
      }
    }
    
    
    if(num.valid.values > 0) {
      #decide the alternative if alternative == "auto"
      if(alternative == "auto") {
        alt <- ifelse(orig.ev < mean(rand.ev, na.rm=TRUE), "less", "greater")
      } else {
        alt <- alternative
      }
      
      #Compute the p-value
      if (alt == "less") {
        pval <- (sum(orig.ev >= rand.ev, na.rm=TRUE) + 1) / (num.valid.values + 1)
      } else { #alt == "greater"
        pval <- (sum(orig.ev <= rand.ev, na.rm=TRUE) + 1) / (num.valid.values + 1)
      }
      #if the original alternative was not the best one, suggest the user to change it
      if(alternative=="greater" & orig.ev<mean(rand.ev,na.rm=TRUE)) message("Alternative is greater and the observed statistic is less than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
      if(alternative=="less" & orig.ev>mean(rand.ev,na.rm=TRUE)) message("Alternative is less and the observed statistic is greater than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
      #Compute the z-score
      if(orig.ev == 0 & all(rand.ev == 0)){ #If everything is 0, warning and "empty" results
        warning(paste0("All permuted values and the original evaluation value are equal to 0. Z-score cannot be computed."))
        pval <- 1
        zscore <- NA
      } else{
        zscore <- round((orig.ev - mean(rand.ev, na.rm=TRUE)) / stats::sd(rand.ev, na.rm=TRUE), 4)
      }
    } else {
      pval <- NA
      zscore <- NA
      alt <- alternative
    }
    
    
    if(!is.na(pval)) {
      if(!is.finite(zscore)){ #if all evaluations are equal, the sd is 0 and the z-score is infinite
        warning(paste0("All permuted values are equal to ", rand.ev[1], ". Z-score is infinite."))
      }  
    }
    
    #Create the permTestResults object
    res<-list(pval=pval, ntimes=ntimes, alternative=alt, observed=orig.ev, permuted=rand.ev, zscore=zscore,
              evaluate.function=evaluate.function[[i]], evaluate.function.name=func.name,
              randomize.function=randomize.function, randomize.function.name=randomize.function.name)
    class(res) <- "permTestResults"  
    results[[func.name]] <- res
  }
  
  class(results) <- "permTestResultsList"
  
  return(results)
  
}
####################################################################################################################
####################################################################################################################
