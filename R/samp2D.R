#' Two Variable Rejection Sampler
#'
#' This function implements two variable rejection sampling for rvs with bounded support and which have bounded pdf.
#'
#' @param fun The pdf that we are sampling from, input as a string. The two variables sampled from must be input using the variables "x" and "y".
#'
#' @param a The lower bound of the variable "x" in the pdf.
#'
#' @param b The upper bound of the variable "x" in the pdf.
#'
#' @param c The lower bound of the variable "y" in the pdf.
#'
#' @param d The upper bound of the variable "y" in the pdf.
#'
#' @param N the number of samples output by the sampling function
#'
#' @return A matrix with N rows and 2 columns containing the samples from the pdf.
#' @export
#'
#' @examples
#'
#' f <- function(z) {
#'  x <- z[1]
#'  y <- z[2]
#'  ifelse(0 <= x & x <= 1 & 0 <= y & y <= 1, x + y, 0)
#' }
#' samps <- samp2D(f, 10000)
#' samps <- data.frame(samps)
#' colnames(samps) <- c("x","y")
#' ggplot(samps, aes(x, y)) +
#'  geom_density_2d()
#'


samp2D <- function(f,N) {
  text <- toString(body(f))
  xminpos <- rep(0,10)
  xmaxpos <- rep(0,10)
  yminpos <- rep(0,10)
  ymaxpos <- rep(0,10)

  # x lower bound
  i <- 1
  pxminpos <- gregexpr('< x ',text)
  for (j in 1:length(pxminpos[[1]])) {
    xminpos[i] <- pxminpos[[1]][j]
    i <- i + 1
  }
  pxminpos <- gregexpr('<= x ',text)
  for (j in 1:length(pxminpos[[1]])) {
    xminpos[i] <- pxminpos[[1]][j]
    i <- i + 1
  }
  pxminpos <- gregexpr('x > ',text)
  for (j in 1:length(pxminpos[[1]])) {
    xminpos[i] <- pxminpos[[1]][j]
    i <- i + 1
  }
  pxminpos <- gregexpr('x >= ',text)
  for (j in 1:length(pxminpos[[1]])) {
    xminpos[i] <- pxminpos[[1]][j]
    i <- i + 1
  }
  xminpos <- xminpos[which(xminpos > 0)]
  xmins <- rep(0,length(xminpos))
  if (length(xmins) != 0) {
    for (j in 1:length(xmins)) {
      xmins[j] <- as.numeric(gsub('[^0-9//-//.]',"",substr(text,xminpos[j]-5,xminpos[j]+5)))
    }
    xmin <- min(xmins[which(!is.na(xmins))])
  }
  else {
    xmin <- NA
  }

  # x upper bound
  i <- 1
  pxmaxpos <- gregexpr('x < ',text)
  for (j in 1:length(pxmaxpos[[1]])) {
    xmaxpos[i] <- pxmaxpos[[1]][j]
    i <- i + 1
  }
  pxmaxpos <- gregexpr('x <= ',text)
  for (j in 1:length(pxmaxpos[[1]])) {
    xmaxpos[i] <- pxmaxpos[[1]][j]
    i <- i + 1
  }
  pxmaxpos <- gregexpr('> x ',text)
  for (j in 1:length(pxmaxpos[[1]])) {
    xmaxpos[i] <- pxmaxpos[[1]][j]
    i <- i + 1
  }
  pxmaxpos <- gregexpr('>= x ',text)
  for (j in 1:length(pxmaxpos[[1]])) {
    xmaxpos[i] <- pxmaxpos[[1]][j]
    i <- i + 1
  }
  xmaxpos <- xmaxpos[which(xmaxpos > 0)]
  xmaxs <- rep(0,length(xmaxpos))
  if (length(xmaxs) != 0) {
    for (j in 1:length(xmaxs)) {
      xmaxs[j] <- as.numeric(gsub('[^0-9//-//.]',"",substr(text,xmaxpos[j]-5,xmaxpos[j]+5)))
    }
    xmax <- max(xmaxs[which(!is.na(xmaxs))])
  }
  else {
    xmax <- NA
  }

  # y lower bound
  i <- 1
  pyminpos <- gregexpr('< y ',text)
  for (j in 1:length(pyminpos[[1]])) {
    yminpos[i] <- pyminpos[[1]][j]
    i <- i + 1
  }
  pyminpos <- gregexpr('<= y ',text)
  for (j in 1:length(pyminpos[[1]])) {
    yminpos[i] <- pyminpos[[1]][j]
    i <- i + 1
  }
  pyminpos <- gregexpr('y > ',text)
  for (j in 1:length(pyminpos[[1]])) {
    yminpos[i] <- pyminpos[[1]][j]
    i <- i + 1
  }
  pyminpos <- gregexpr('y >= ',text)
  for (j in 1:length(pyminpos[[1]])) {
    yminpos[i] <- pyminpos[[1]][j]
    i <- i + 1
  }
  yminpos <- yminpos[which(yminpos > 0)]
  ymins <- rep(0,length(yminpos))
  if (length(ymins) != 0) {
    for (j in 1:length(ymins)) {
      ymins[j] <- as.numeric(gsub('[^0-9//-//.]',"",substr(text,yminpos[j]-5,yminpos[j]+5)))
    }
    ymin <- min(ymins[which(!is.na(ymins))])
  }
  else {
    ymin <- NA
  }

  # y upper bound
  i <- 1
  pymaxpos <- gregexpr('y < ',text)
  for (j in 1:length(pymaxpos[[1]])) {
    ymaxpos[i] <- pymaxpos[[1]][j]
    i <- i + 1
  }
  pymaxpos <- gregexpr('y <= ',text)
  for (j in 1:length(pymaxpos[[1]])) {
    ymaxpos[i] <- pymaxpos[[1]][j]
    i <- i + 1
  }
  pymaxpos <- gregexpr('> y ',text)
  for (j in 1:length(pymaxpos[[1]])) {
    ymaxpos[i] <- pymaxpos[[1]][j]
    i <- i + 1
  }
  pymaxpos <- gregexpr('>= y ',text)
  for (j in 1:length(pymaxpos[[1]])) {
    ymaxpos[i] <- pymaxpos[[1]][j]
    i <- i + 1
  }
  ymaxpos <- ymaxpos[which(ymaxpos > 0)]
  ymaxs <- rep(0,length(ymaxpos))
  if (length(ymaxs) != 0){
    for (j in 1:length(ymaxs)) {
      ymaxs[j] <- as.numeric(gsub('[^0-9//-//.]',"",substr(text,ymaxpos[j]-5,ymaxpos[j]+5)))
    }
    ymax <- max(ymaxs[which(!is.na(ymaxs))])
  }
  else {
    ymax <- NA
  }

  samples <- matrix(rep(0,2*N), nrow = N, ncol = 2) # creating a matrix to store the samples in
  if (!is.na(xmin) & !is.na(xmax) & !is.na(ymin) & !is.na(ymax)) {
    maxf <- optim(c(xmax/4,ymax/4),f,control = list(fnscale = -1)) # finds the maximum of the pdf
    maxf <- maxf$value + .1
    i <- 0
    while( i < N) { # while loop that runs unil N accepted samples are found.
      psx <- runif(1,xmin,xmax) # creating a potential x value to be tested
      psy <- runif(1,ymin,ymax) # creating a potential y value to be tested
      testsamp <- runif(1,0, maxf) # creating a test sample using maxf$value to find the max of the pdf and adding .1 to account for errors in tolerance.
      if (testsamp < f(c(psx,psy))) { # tests if the potential x and y values should be rejected. If kept, the x value is stored in the first column and the y value in the second of the samples matrix.
        samples[i+1,1] <- psx
        samples[i+1,2] <- psy
        i <- i +1
      }
    }
  }
  else {
    k <- 1
    maxfvalues <- rep(0,6561)
    maxs <- matrix(rep(0,6561*2),ncol = 2, byrow = TRUE)
    for (j in seq(-10,10,.25)) {
      for (i in seq(-10,10,.25)) {
        value <- optim(c(j,i),f,control = list(fnscale = -1))
        maxfvalues[k] <- value$value
        maxs[k,] <- value$par
        k <- k + 1
      }
    }
    maxf <- max(maxfvalues)
    mean <- maxs[which(maxfvalues == maxf),]
    if (sum(dim(mean))>2) {
      mean <- mean[1,]
    }
    check <- rep(0,6561)
    while (sum(check) != 6561) {
      check <- rep(0,6561)
      k <- 1
      for (j in seq(-10,10,.25)) {
        for (i in seq(-10,10,.25)) {
         if (maxf^2*dnorm(i,mean[1],5)*dnorm(j,mean[2],5) >= f(c(i,j))) {
           check[k] <- 1
         }
          k <- k + 1
        }
      }
      maxf <- maxf + 1
    }
    maxf <- maxf - 1
    i <- 0
    while (i < N) {
      psx <- rnorm(1,mean[1],5)
      psy <- rnorm(1,mean[2],5)
      testsamp <- runif(1,0,maxf^2*dnorm(psx,mean[1],5)*dnorm(psy,mean[2],5))
      if (testsamp < f(c(psx,psy))) {
        samples[i+1,1] <- psx
        samples[i+1,2] <- psy
        i <- i + 1
      }
    }
  }
  samples # outputs the samples matrix.
}
