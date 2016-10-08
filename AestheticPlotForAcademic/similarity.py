simiheat <- function(labs, ...){
  library('ggplot2')
  dots <- list(...)
  cols = length(dots)
  rInterset <- c()
  for (outterIndex in seq_along(dots)){
    for (innerIndex in seq_along(dots)){
      rInterset <- c(rInterset, length(intersect(dots[[outterIndex]], dots[[innerIndex]]))/length(dots[[innerIndex]]))
    }
  }
  y <- rep(1:cols, cols)
  x <- rep(1:cols, each = cols)
  df <- as.data.frame(cbind(x, y, rInterset))
  draftPlot <- ggplot(df, aes(x, y))
  draftPlot+geom_tile(aes(fill=rInterset))
    +scale_y_continuous(expand=c(0, 0), breaks = c(1:cols), labels=labs)
    +scale_x_continuous(expand=c(0, 0), breaks = c(1:cols), labels=labs)
    +scale_fill_gradientn(limits=c(0, 1), colours=c("white", "black"))
    +theme_linedraw()+theme(axis.text.x=element_text(angle = 90))
    +guides(fill = guide_legend(title = NULL))
  draftPlot
}
