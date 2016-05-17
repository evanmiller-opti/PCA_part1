# build example where even and odd variables are bringing in noisy images
# of two different signals.
mkData <- function(n) {
  for(group in 1:10) {
    # y is the sum of two effects yA and yB
    yA <- rnorm(n)
    yB <- rnorm(n)
    if(group==1) {
      d <- data.frame(y=yA+yB+rnorm(n))
      code <- 'x'
    } else {
      code <- paste0('noise',group-1)
    }
    yS <- list(yA,yB)
    # these variables are correlated with y in group 1,
    # but only to each other (and not y) in other groups
    for(i in 1:5) {
      vi <- yS[[1+(i%%2)]] + rnorm(nrow(d))
      d[[paste(code,formatC(i,width=2,flag=0),sep='.')]] <- ncol(d)*vi
    }
  }
  d
}

haveWVPlots <- require('WVPlots',quietly=TRUE) # https://github.com/WinVector/WVPlots

barbell_plot = function(frame, xvar, ymin, ymax, colorvar=NULL) {
  if(is.null(colorvar)) {
    gplot = ggplot(frame, aes_string(x=xvar))
  } else {
    gplot = ggplot(frame, aes_string(x=xvar, color=colorvar))
  }
  
  gplot + geom_point(aes_string(y=ymin)) + 
    geom_point(aes_string(y=ymax)) +
    geom_linerange(aes_string(ymin=ymin, ymax=ymax))
}

dotplot_identity = function(frame, xvar, yvar, colorvar=NULL) {
  if(is.null(colorvar)) {
    gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar))
  } else {
    gplot = ggplot(frame, aes_string(x=xvar, y=yvar, ymax=yvar, color=colorvar))
  }
  
  gplot + geom_point() + geom_linerange(ymin=0)
}

extractProjection <- function(ndim,princ) {
  # pull off the rotation.  
  proj <- princ$rotation[,1:ndim] 
  # sign was arbitrary, so flip in convenient form
  for(i in seq_len(ndim)) {
    si <- sign(mean(proj[,i]))
    if(si!=0) {
      proj[,i] <- proj[,i]*si
    }
  }
  proj
}

rsq <- function(pred,y) {
  1-sum((y-pred)^2)/sum((y-mean(y))^2)
}