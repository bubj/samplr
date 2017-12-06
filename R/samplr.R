#' Single and Two Variable Rejection Sampler
#'
#' This function implements both one and two variable rejection sampling for probability density functions which are bounded, but the support does not need to be bounded.
#'
#' @param f The pdf that we are sampling from, input as a predefined function. For one variable, the input variable must be declared as \code{x}. For two variables, the parameter of the function must be a vector containing the variables \code{x} and \code{y}.
#'
#' @param N The number of samples returned by the rejection sampler.
#'
#' @param twod Specification if the input pdf is two dimensional if \code{TRUE}. Default is \code{FALSE}, meaning that the input function is one dimensional.
#'
#' @return If using a single variable, the output is a vector containing the samples. If two variables, the output is a matrix with N rows and 2 columns containing the samples.
#' @export
#'
#' @examples
#'
#' Single variable
#'
#' f <- function(x) {
#'   ifelse(0 < x & x < 1, 2*x, 0)
#' }
#' hist(samplr(f,10000))
#'
#' g <- function(x) {
#'   ifelse(0 < x & x < 2, 1/2*x, 0)
#' }
#' hist(samplr(g,10000))
#'
#' h <- function(x) {
#'   ifelse(0 < x & x < 6.2832, 1/2/pi*(sin(x) + 1), 0)
#' }
#' hist(samplr(h,10000))
#'
#' u <- function(x) {
#'   ifelse(x < 0, 1/2*exp(x), 1/2*exp(-x))
#' }
#' hist(samplr(u,10000))
#'
#' Two dimensional (requires \code{ggplot2})
#'
#' v <- function(z) {
#'   x <- z[1]
#'   y <- z[2]
#'   ifelse(0 <= x & x <= 1 & 0 <= y & y <= 1, x + y, 0)
#' }
#' samps <- data.frame(samplr(v, 10000, twod = TRUE))
#' colnames(samps) <- c("x","y")
#' ggplot(samps, aes(x, y)) +
#'   geom_density_2d()
#'
#' w <- function(z) {
#'   x <- z[1]
#'   y <- z[2]
#'   ifelse(0 <= x & x <= 1 & 0 <= y & y <= 1 & 0 <= x + y & x + y <= 1, 24*x*y,0)
#' }
#' samps <- data.frame(samplr(w, 10000, twod = TRUE))
#' colnames(samps) <- c("x","y")
#' ggplot(samps, aes(x, y)) +
#'   geom_density_2d()
#'
#' l <- function(z) {
#'   x <- z[1]
#'   y <- z[2]
#'   ifelse(0 < x & 0 < y, exp(-x)*exp(-y), 0)
#' }
#' samps <- data.frame(samplr(l, 10000, twod = TRUE))
#' colnames(samps) <- c("x","y")
#' ggplot(samps, aes(x, y)) +
#'   geom_density_2d()
#'


samplr <- function(f, N, twod = FALSE) {
  if (twod == FALSE) {
    text <- toString(body(f)) # convert body of the function to a string
    leftboundpos <- rep(0,10)
    rightboundpos <- rep(0,10)

    # find the left bound of the function by searching for where the variable x is greater than some value
    i <- 1
    pleftpos <- gregexpr('< x',text)
    for (j in 1:length(pleftpos[[1]])) {
      leftboundpos[i] <- pleftpos[[1]][j]
      i <- i + 1
    }
    pleftpos <- gregexpr('<= x',text)
    for (j in 1:length(pleftpos[[1]])) {
      leftboundpos[i] <- pleftpos[[1]][j]
      i <- i + 1
    }
    pleftpos <- gregexpr('x > ',text)
    for (j in 1:length(pleftpos[[1]])) {
      leftboundpos[i] <- pleftpos[[1]][j]
      i <- i + 1
    }
    pleftpos <- gregexpr('x >= ',text)
    for (j in 1:length(pleftpos[[1]])) {
      leftboundpos[i] <- pleftpos[[1]][j]
      i <- i + 1
    }
    leftboundpos <- leftboundpos[which(leftboundpos > 0)]
    leftbounds <- rep(0,length(leftboundpos))
    if (length(leftboundpos) != 0) {
      for (j in 1:length(leftboundpos)) {
        leftbounds[j] <- as.numeric(gsub("[^0-9\\-\\.]","",(substr(text,leftboundpos[j]-6,leftboundpos[j]+6))))
      }
      leftbound <- min(leftbounds)
    }
    else {
      leftbound <- NA
    }

    #find the right bound of the function by searching where the variable x is less than a value
    i <- 1
    prightpos <- gregexpr('> x',text)
    for (j in 1:length(prightpos[[1]])) {
      rightboundpos[i] <- prightpos[[1]][j]
      i <- i + 1
    }
    prightpos <- gregexpr('>= x',text)
    for (j in 1:length(prightpos[[1]])) {
      rightboundpos[i] <- prightpos[[1]][j]
      i <- i + 1
    }
    prightpos <- gregexpr('x < ',text)
    for (j in 1:length(prightpos[[1]])) {
      rightboundpos[i] <- prightpos[[1]][j]
      i <- i + 1
    }
    prightpos <- gregexpr('x <= ',text)
    for (j in 1:length(prightpos[[1]])) {
      rightboundpos[i] <- prightpos[[1]][j]
      i <- i + 1
    }
    rightboundpos <- rightboundpos[which(rightboundpos > 0)]
    rightbounds <- rep(0,length(rightboundpos))
    if (length(rightboundpos) != 0) {
      for (j in 1:length(rightboundpos)) {
        rightbounds[j] <- as.numeric(gsub("[^0-9\\-\\.]","",(substr(text,rightboundpos[j]-6,rightboundpos[j]+6))))
      }
      rightbound <- max(rightbounds)
    }
    else {
      rightbound <- NA
    }

    if (!is.na(rightbound) & !is.na(leftbound)) {
      if ( length(rightbounds) == length(leftbounds) ) {
        if ( length(rightbounds) == 1 ) {
          if ( rightbounds[1] == leftbounds[1] ) {
            rightbound <- NA
            leftbound <- NA
          }
          else {
            rightbound <- max(rightbounds)
            leftbound <- min(leftbounds)
          }
        }
        else {
          rightbound <- max(rightbounds)
          leftbound <- min(leftbounds)
        }
      }
      else {
        if (length(rightbounds) < length(leftbounds)) {
          rightbound <- NA
          leftbound <- min(leftbounds)
        }
        else {
          rightbound <- max(rightbounds)
          leftbound <- NA
        }
      }
    }

    samples <- rep(0,N) # creates a vector for storing the sample values
    if (!is.na(rightbound) & !is.na(leftbound)) {
      maxf <- optimize(f,c(leftbound,rightbound),maximum = TRUE)
      maxf <- maxf$objective + .1
      i <- 0
      while ( i < N) { # while loop that runs until N samples are selected
        potsamp <- runif(1,leftbound,rightbound) # creating a potential sample
        testsamp <- runif(1, 0, maxf) # creating a test sample. Uses maxf$objective to get the max of the function, and adds .1 to account for errors in tolerance
        if ( testsamp < f(potsamp) ) { # tests if the potential sample should be rejected. If it is kept, it is kept in the samples vector
          samples[i+1] = potsamp
          i = i + 1
        }
      }
    }

    else {
      maxfvalues <- rep(0,200)
      maxes <- rep(0,200)
      for (i in -100:99) {
        value <- optimize(f,c(i-.1,i+1.1),maximum = TRUE)
        maxfvalues[i+101] <- value$objective
        maxes[i+101] <- value$maximum
      }
      maxf <- max(maxfvalues)
      mean <- max(maxes[which(maxfvalues == maxf)])

      check <- rep(0,101)
      while (sum(check) != 101) {
        check <- rep(0,101)
        for (i in -50:50) {
          if ( suppressWarnings(maxf*dt(i,1,mean)) >= f(i) ) {
            check[i+51] <- 1
          }
        }
        maxf <- maxf + 1
      }
      maxf <- maxf - 1

      i <- 0
      while (i < N){
        potsamp <- rt(1,1,mean)
        testsamp <- runif(1, 0, suppressWarnings(maxf*dt(potsamp,1,mean)))
        if (testsamp < f(potsamp)) {
          samples[i+1] <- potsamp
          i <- i + 1
        }
      }
    }
    samples
  }
  else {
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
}
