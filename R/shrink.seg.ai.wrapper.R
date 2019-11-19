#' unite neighbouring segments if possible
#'
#' @param seg segmentation data
shrink.seg.ai.wrapper<-function(seg){
  new.seg <- seg[1,]
  #For each of  the samples
  for(j in unique(seg[,1])){
    sample.seg <- seg[seg[,1] %in% j,]
    new.sample.seg <- seg[1,]
    #For each of the chromosomes
    for(i in unique(sample.seg[,2])){
      sample.chrom.seg <- sample.seg[sample.seg[,2] %in% i,,drop=F]
      #Just make shrink for each chromosomes, with more than two segments.
      if(nrow(sample.chrom.seg) > 1){
        sample.chrom.seg <- shrink.seg.ai(sample.chrom.seg)
      }
      #Concatenate the rows, of each chromosome
      new.sample.seg <- rbind(new.sample.seg, sample.chrom.seg)
    }
    #Concatenate the rows of each sample
    new.seg <- rbind(new.seg, new.sample.seg[-1,])
  }
  seg <- new.seg[-1,]
  return(seg)
}
