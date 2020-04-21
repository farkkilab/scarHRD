#' Preprocessing for further analysis
#'
#' @param seg segmentation data
#' @return preprocessed data
#This function will select the segments bigger in A tmp[,8] > tmp[,7] that do not spam all the chromosome !all(seg[,8] <= seg[,7]
#After reducing it get the number of events using shrink.seg.ai.wrapper
#Important to note, is that columns seg[,7] and seg[,8] are interchanged

preprocess.hrd<-function(seg){
  #Will ignore chromosomes X,Y
  seg <- seg[!seg[,2] %in% c(paste('chr',c('X','Y','x','y',23,24),sep=''),c('X','Y','x','y',23,24)),]
  seg[,1] <- as.character(seg[,1])

  if(! all(seg[,8] <= seg[,7]) ){
    tmp <- seg
    seg[tmp[,8] > tmp[,7],7]  <- tmp[tmp[,8] > tmp[,7],8]
    seg[tmp[,8] > tmp[,7],8]  <- tmp[tmp[,8] > tmp[,7],7]
  }
  seg <- shrink.seg.ai.wrapper(seg)
  return(seg)

}
