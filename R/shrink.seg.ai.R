#' unite neighbouring segments if possible
#'
#' @param chr.seg segmentation data
#' @return gc normalized copy-number
#Example of shrink:
#chr1 1 1000 2 3 2
#chr1 1002 2002 2 3 2
#chr1 2005 3003 2 3 2
#result of shrink:
#chr1 1 3003 2 9 2

shrink.seg.ai<-function(chr.seg){
  new.chr <- chr.seg
  if(nrow(chr.seg) > 1){
    new.chr <- chr.seg
    seg.class <- c(1)
    for(j in 2:nrow(new.chr)){
      #Check if the cnv values for A and B allele in the previous segment is equal in the current segment
      seg_test <- new.chr[(j-1),7] == new.chr[j,7] & new.chr[(j-1),8] == new.chr[j,8]
      if(seg_test){
	  #If equal, then save the rownumber of the previous segment
          seg.class <- c(seg.class, seg.class[j-1])
        }else{	
    seg.class <- c(seg.class, seg.class[j-1]+1)
  }
    }
    #For each of the unique rownumbers (consecutive segments, initial row), select the maximun End_position (final position of the segment)
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
