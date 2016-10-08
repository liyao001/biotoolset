simiheat <- function(..., saveTo, tags, startColor = 'white', endColor = 'black', ang = 90, format = 'pdf', jaccard = TRUE, width = 8, height = 8){
  library('ggplot2')
  dots <- list(...)
  cols = length(dots)
  bs <- c(1:cols)
  rInterset <- c()
  
  if (ang == 90){
    hj <- 1
    vj <- .5
  }else if(ang == 0){
    hj <- NULL
    vj <- NULL
  }else{
    hj <- 1
    vj <- 1
  }
  
  #calculate similarity
  if (jaccard == TRUE){
    for (outterIndex in seq_along(dots)){
      for (innerIndex in seq_along(dots)){
        #Calculate Jaccard index
        rInterset <- c(rInterset, length(intersect(dots[[outterIndex]], dots[[innerIndex]]))/length(union(dots[[outterIndex]], dots[[innerIndex]])))
      }
    }
  }else{
    for (outterIndex in seq_along(dots)){
      for (innerIndex in seq_along(dots)){
        rInterset <- c(rInterset, length(intersect(dots[[outterIndex]], dots[[innerIndex]]))/length(dots[[innerIndex]]))
      }
    }
  }
  
  y <- rep(1:cols, cols)
  x <- rep(1:cols, each = cols)
  df <- as.data.frame(cbind(x, y, rInterset))
  
  #create heatmap
  draftPlot <- ggplot(df, aes(x, y))
  draftPlot <- draftPlot+geom_tile(aes(fill=rInterset))+labs(x=NULL, y=NULL)+scale_y_continuous(expand=c(0, 0), breaks = bs, labels=tags)+scale_x_continuous(expand=c(0, 0), breaks = bs, labels=tags)+scale_fill_gradientn(limits=c(0, 1), colours=c(startColor, endColor))+theme_linedraw()+theme(axis.text.x=element_text(angle = ang, hjust = hj, vjust = vj))+guides(fill = guide_legend(title = NULL))
  #save plot
  ggsave(saveTo, draftPlot, device = format, width = width, height = height)
}