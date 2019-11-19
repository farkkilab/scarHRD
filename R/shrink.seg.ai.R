#' unite neighbouring segments if possible
#'
#' @param chr.seg segmentation data
#' @return gc normalized copy-number
shrink.seg.ai<-function(chr.seg){
  if(nrow(chr.seg) > 1){
    new.chr <- chr.seg
    seg.class <- c(1)
    #Check if the cnv values for A and B allele in the previous segment is equal in the current segment
    for(j in 2:nrow(new.chr)){
      seg_test <- new.chr[(j-1),7] == new.chr[j,7] & new.chr[(j-1),8] == new.chr[j,8]
      if(seg_test){
          seg.class <- c(seg.class, seg.class[j-1])
        }else{
    seg.class <- c(seg.class, seg.class[j-1]+1)
  }
    }
    for(j in unique(seg.class)){
      new.chr[seg.class %in% j,4] <- max(new.chr[seg.class %in% j,4])
      new.chr[seg.class %in% j,5] <- sum(new.chr[seg.class %in% j,5])
    }
    new.chr<- new.chr[!duplicated(seg.class),]
  }
  if(nrow(chr.seg) == 1){
    new.chr <- chr.seg
  }
  return(new.chr)
}
