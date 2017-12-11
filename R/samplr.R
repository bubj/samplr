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
#' hist(samplr(10000,f))
#'
#' g <- function(x) {
#'   ifelse(0 < x & x < 2, 1/2*x, 0)
#' }
#' hist(samplr(10000,g))
#'
#' h <- function(x) {
#'   ifelse(0 < x & x < 6.2832, 1/2/pi*(sin(x) + 1), 0)
#' }
#' hist(samplr(10000,h))
#'
#' u <- function(x) {
#'   ifelse(x < 0, 1/2*exp(x), 1/2*exp(-x))
#' }
#' hist(samplr(10000,u))
#'
#' Two dimensional (requires ggplot2)
#'
#' v <- function(z) {
#'   x <- z[1]
#'   y <- z[2]
#'   ifelse(0 <= x & x <= 1 & 0 <= y & y <= 1, x + y, 0)
#' }
#' samps <- data.frame(samplr(10000, v, twod = TRUE))
#' colnames(samps) <- c("x","y")
#' ggplot(samps, aes(x, y)) +
#'   geom_density_2d()
#'
#' w <- function(z) {
#'   x <- z[1]
#'   y <- z[2]
#'   ifelse(0 <= x & x <= 1 & 0 <= y & y <= 1 & 0 <= x + y & x + y <= 1, 24*x*y,0)
#' }
#' samps <- data.frame(samplr(10000, w, twod = TRUE))
#' colnames(samps) <- c("x","y")
#' ggplot(samps, aes(x, y)) +
#'   geom_density_2d()
#'
#' l <- function(z) {
#'   x <- z[1]
#'   y <- z[2]
#'   ifelse(0 < x & 0 < y, exp(-x)*exp(-y), 0)
#' }
#' samps <- data.frame(samplr(10000, l, twod = TRUE))
#' colnames(samps) <- c("x","y")
#' ggplot(samps, aes(x, y)) +
#'   geom_density_2d()
#'


samplr <- function(N, f, twod = FALSE) {
  if (twod == FALSE) {
    text <- toString(body(f)) # convert body of the function to a string
    # create vectors to store positions of the string where bounds are defined
    leftboundpos <- rep(0,10)
    rightboundpos <- rep(0,10)

    # find the left bound of the function by searching for where the variable x is greater than some value
    i <- 1
    # search the body of the function for '< x', and store its position in the vector.
    pleftpos <- gregexpr('< x',text)
    for (j in 1:length(pleftpos[[1]])) {
      leftboundpos[i] <- pleftpos[[1]][j]
      i <- i + 1
    }
    # search the body of the function for '<= x', and store its position in the vector.
    pleftpos <- gregexpr('<= x',text)
    for (j in 1:length(pleftpos[[1]])) {
      leftboundpos[i] <- pleftpos[[1]][j]
      i <- i + 1
    }
    # search the body of the function for 'x >', and store its position in the vector.
    pleftpos <- gregexpr('x > ',text)
    for (j in 1:length(pleftpos[[1]])) {
      leftboundpos[i] <- pleftpos[[1]][j]
      i <- i + 1
    }
    # search the body of the function for 'x >= ', and store its position in the vector.
    pleftpos <- gregexpr('x >= ',text)
    for (j in 1:length(pleftpos[[1]])) {
      leftboundpos[i] <- pleftpos[[1]][j]
      i <- i + 1
    }
    leftboundpos <- leftboundpos[which(leftboundpos > 0)] # keep only the positions found that are not less than 0
    leftbounds <- rep(0,length(leftboundpos)) # create a vector for possible lower bounds
    if (length(leftboundpos) != 0) { # check that there are some positions that have been found that are greater than 0
      for (j in 1:length(leftboundpos)) { # iterate through all of the positions, keeping only numbers in the vicinity of the position found. These numbers are the lower bounds.
        leftbounds[j] <- as.numeric(gsub("[^0-9\\-\\.]","",(substr(text,leftboundpos[j]-6,leftboundpos[j]+6))))
      }
      leftbound <- min(leftbounds) # the lower bound is the min of the values found through the iteration.
    }
    else { # if there were no positions found for the lower bounds, then the lower bound is not defined.
      leftbound <- NA
    }

    #find the right bound of the function by searching where the variable x is less than a value
    i <- 1
    # search the body of the function for '> x', and store its position in the vector.
    prightpos <- gregexpr('> x',text)
    for (j in 1:length(prightpos[[1]])) {
      rightboundpos[i] <- prightpos[[1]][j]
      i <- i + 1
    }
    # search the body of the function for '>= x', and store its position in the vector.
    prightpos <- gregexpr('>= x',text)
    for (j in 1:length(prightpos[[1]])) {
      rightboundpos[i] <- prightpos[[1]][j]
      i <- i + 1
    }
    # search the body of the function for 'x < ', and store its position in the vector.
    prightpos <- gregexpr('x < ',text)
    for (j in 1:length(prightpos[[1]])) {
      rightboundpos[i] <- prightpos[[1]][j]
      i <- i + 1
    }
    # search the body of the function for 'x <= ', and store its position in the vector.
    prightpos <- gregexpr('x <= ',text)
    for (j in 1:length(prightpos[[1]])) {
      rightboundpos[i] <- prightpos[[1]][j]
      i <- i + 1
    }
    rightboundpos <- rightboundpos[which(rightboundpos > 0)] # keep only positions that are greater than 0
    rightbounds <- rep(0,length(rightboundpos)) # create a vector for possible upper bounds
    if (length(rightboundpos) != 0) { # check that there were positive positions found
      for (j in 1:length(rightboundpos)) {
        rightbounds[j] <- as.numeric(gsub("[^0-9\\-\\.]","",(substr(text,rightboundpos[j]-6,rightboundpos[j]+6)))) # iterate through all the positions, keeping only numbers in the string in the vicinity of each position.
      }
      rightbound <- max(rightbounds) # the upper bound is the max of the bounds found through iteration.
    }
    else { # if there are no positions greater than 0, then the upper bound is undefined.
      rightbound <- NA
    }

    if (!is.na(rightbound) & !is.na(leftbound)) {
      if ( length(rightbounds) == length(leftbounds) ) {
        if ( length(rightbounds) == 1 ) {
          if ( rightbounds[1] == leftbounds[1] ) { # if there was only one upper and lower bound found, and they are equal to each other, then the bound of the support is undefined
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
        if (length(rightbounds) < length(leftbounds)) { # if the total number of upper bounds found is less than the number of lower bounds found, then the upper bound is undefined
          rightbound <- NA
          leftbound <- min(leftbounds)
        }
        else { # if the total number of lower bounds found is less than the number of upper bounds, then the lower bound is undefined
          rightbound <- max(rightbounds)
          leftbound <- NA
        }
      }
    }

    samples <- rep(0,N) # creates a vector for storing the sample values
    if (!is.na(rightbound) & !is.na(leftbound)) { # if the bounds are defined
      maxf <- optimize(f,c(leftbound,rightbound),maximum = TRUE) # find the maximum of the function
      maxf <- maxf$objective + .1 # add .1 to account for tolerance
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

    else { # if one of the bounds is undefined
      maxfvalues <- rep(0,200) # create a vector to store maximum values of the function
      maxes <- rep(0,200) # create a vector to store locations of maximum values
      for (i in -100:99) { # try to find maxmiums in this range.
        value <- optimize(f,c(i-.1,i+1.1),maximum = TRUE)
        maxfvalues[i+101] <- value$objective
        maxes[i+101] <- value$maximum
      }
      maxf <- max(maxfvalues) # the maximum value of the pdf is the max of the values found.
      mean <- max(maxes[which(maxfvalues == maxf)]) # the location of the max is stored, which is used as the non-centrality parameter in the bounding t-disribution

      checks <- c(-10000000,-1000000,-100000,-10000,-1000,-500,-250,-100:100,250,500,1000,10000,100000,1000000,10000000) # these x-values are checked to make sure the bounding t-distribution is above the pdf
      check <- rep(0,length(checks))
      while (sum(check) != length(check)) { # while there are some values in the range set where the t distribution is below the pdf, iterate using increasing maxf values
        check <- rep(0,length(checks))
        i <- 1
        for (i in 1:length(check)) {
          if ( suppressWarnings(maxf*dt(checks[i],1,mean)) >= f(checks[i])) { # check to see that the t distribution is above the pdf
            check[i] <- 1
          }
        }
        maxf <- maxf + 1
      }

      i <- 0
      while (i < N){
        potsamp <- rt(1,1,mean) # create a potential sample
        testsamp <- runif(1, 0, suppressWarnings(maxf*dt(potsamp,1,mean))) # create a test sample
        if (testsamp < f(potsamp)) { # check that the test sample is below the pdf
          samples[i+1] <- potsamp
          i <- i + 1
        }
      }
    }
    samples # output the samples
  }
  else { # if a 2D pdf is specified
    text <- toString(body(f)) # create a string from the body of the function, and create vectors to store mins and maxes of the support
    xminpos <- rep(0,10)
    xmaxpos <- rep(0,10)
    yminpos <- rep(0,10)
    ymaxpos <- rep(0,10)

    # x lower bound
    i <- 1
    # search the body of the function for '< x ', and store its position in the vector.
    pxminpos <- gregexpr('< x ',text)
    for (j in 1:length(pxminpos[[1]])) {
      xminpos[i] <- pxminpos[[1]][j]
      i <- i + 1
    }
    # search the body of the function for '<= x ', and store its position in the vector.
    pxminpos <- gregexpr('<= x ',text)
    for (j in 1:length(pxminpos[[1]])) {
      xminpos[i] <- pxminpos[[1]][j]
      i <- i + 1
    }
    # search the body of the function for 'x > ', and store its position in the vector.
    pxminpos <- gregexpr('x > ',text)
    for (j in 1:length(pxminpos[[1]])) {
      xminpos[i] <- pxminpos[[1]][j]
      i <- i + 1
    }
    # search the body of the function for 'x >= ', and store its position in the vector.
    pxminpos <- gregexpr('x >= ',text)
    for (j in 1:length(pxminpos[[1]])) {
      xminpos[i] <- pxminpos[[1]][j]
      i <- i + 1
    }
    xminpos <- xminpos[which(xminpos > 0)] # remove all non-negative positions
    xmins <- rep(0,length(xminpos)) # create a vector to store values of the bound
    if (length(xmins) != 0) { # if there are positions that are greater than 0
      for (j in 1:length(xmins)) {
        xmins[j] <- as.numeric(gsub('[^0-9//-//.]',"",substr(text,xminpos[j]-5,xminpos[j]+5))) # remove everything is the string body around each position besides numbers
      }
      xmin <- min(xmins[which(!is.na(xmins))]) # the minimum is the min of the values found
    }
    else {
      xmin <- NA # if no bounds are found, then the bound is undefined
    }

    # x upper bound in the same fashion as x lower bound
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

    # y lower bound in the same fashion as x lower bound
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

    # y upper bound in the same fashion as x lower bound
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
    #print(xmins)
    #print(xmaxs)
    #print(ymins)
    #print(ymaxs)

    # if (!is.na(xmax) & !is.na(xmin)) {
    #   if ( length(xmaxs) == length(xmins) ) {
    #     if ( length(xmaxs) == 1 ) {
    #       if ( xmaxs[1] == xmins[1] ) {
    #         xmax <- NA
    #         xmin <- NA
    #       }
    #       else {
    #         xmax <- max(xmaxs)
    #         xmin <- min(xmins)
    #       }
    #     }
    #     else {
    #       xmax <- max(xmaxs)
    #       xmin <- min(xmins)
    #     }
    #   }
    #   else {
    #     if (length(xmaxs) < length(xmins)) {
    #       xmax <- NA
    #       xmin <- min(xmins)
    #     }
    #     else {
    #       xmax <- max(xmaxs)
    #       xmin <- NA
    #     }
    #   }
    # }
    #
    # if (!is.na(ymax) & !is.na(ymin)) {
    #   if ( length(ymaxs) == length(ymins) ) {
    #     if ( length(ymaxs) == 1 ) {
    #       if ( ymaxs[1] == ymins[1] ) {
    #         ymax <- NA
    #         ymin <- NA
    #       }
    #       else {
    #         ymax <- max(ymaxs)
    #         ymin <- min(ymins)
    #       }
    #     }
    #     else {
    #       ymax <- max(ymaxs)
    #       ymin <- min(ymins)
    #     }
    #   }
    #   else {
    #     if (length(ymaxs) < length(ymins)) {
    #       ymax <- NA
    #       ymin <- min(ymins)
    #     }
    #     else {
    #       ymax <- max(ymaxs)
    #       ymin <- NA
    #     }
    #   }
    # }
    #print(c(xmin,xmax,ymin,ymax))
    # if (length(xmins) != length(xmaxs)) {
    #   xmin <- NA
    #   xmax <- NA
    # }
    # if (length(ymins) != length(ymaxs)) {
    #   ymin <- NA
    #   ymax <- NA
    # }
    if (!is.na(xmin) & !is.na(xmax)) { # if the lower and upper bounds are defined as the same value, the support is unbounded
      if (xmin == xmax) {
        xmin <- NA
        xmax <- NA
      }
    }
    if (!is.na(ymin) & !is.na(ymax)) {
      if (ymin == ymax) {
        ymin <- NA
        ymax <- NA
      }
    }
    #print(c(xmin,xmax,ymin,ymax))

    samples <- matrix(rep(0,2*N), nrow = N, ncol = 2) # creating a matrix to store the samples in
    if (!is.na(xmin) & !is.na(xmax) & !is.na(ymin) & !is.na(ymax)) { # if the support is bounded
      maxf <- optim(c(xmax/4,ymax/4),f,control = list(fnscale = -1)) # finds the maximum of the pdf
      maxf <- maxf$value + .1 # add .1 to account for tolerance errors
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
    else { # if the support is not bounded
      side <- seq(-20,20,.5) # create a vector of number from -20 to 20 by steps of .5
      maxpoints <- matrix(rep(0,2*(length(side)^2)),ncol = 2) # create a matrix of dimenions of the cartesian product of side
      k <- 1
      for(i in 1:length(side)) { # store the cartesian product of the vector side in the matrix
        for(j in 1:length(side)) {
          maxpoints[k,1] <- side[i]
          maxpoints[k,2] <- side[j]
          k <- k + 1
        }
      }
      maxfvalues <- seq(0,length(side)^2-1) # create a vector to store maximum values
      maxes <- matrix(rep(0,2*(length(side)^2)),ncol = 2) # create a matrix to store points of the maximums
      for (i in 1:length(maxfvalues)) { # for each number in the vector of maxfvalues
        value <- optim(c(maxpoints[i,1],maxpoints[i,2]),f,control = list(fnscale = -1)) # find the max
        maxfvalues[i] <- value$value # store the max
        maxes[i,] <- value$par # store the points of the max
      }
      maxf <- max(maxfvalues) # find the largest max
      mean <- maxes[which(maxfvalues == maxf),] # use the points that correspond to the largest max as non-centrality parameters in the r distribution
      #print(maxf)
      #print(mean)
      #print(dim(mean))
      if (!is.null(dim(mean))) { # if there are more than one set of points that correspond to the max value, pick the first
        mean <- mean[1,]
      }

      side <- c(-10000000,-1000000,-100000,-10000,-1000,-500,-250,seq(-100,100,5),250,500,1000,10000,100000,1000000,10000000) # create a vector for points that need to be checked

      # create a matrix and store the cartesian product of the side vector in the matrix
      checks <- matrix(rep(0,2*(length(side)^2)),ncol = 2)
      k <- 1
      for (i in 1:length(side)) {
        for (j in 1:length(side)) {
          checks[k,1] <- side[i]
          checks[k,2] <- side[j]
          k <- k + 1
        }
      }
      check <- rep(0,dim(checks)[1])
      while (sum(check) != length(check)) { # while there are still points in the t distribution that are below the pdf, iterate with increasing maxf
        check <- rep(0,dim(checks)[1])
        for (i in 1:length(check)) {
          if (suppressWarnings(maxf*dt(checks[i,1],1,mean[1])*dt(checks[i,2],1,mean[2])) >= f(c(checks[i,1],checks[i,2]))) # check that the product of t distrubtions is above the pdf
            check[i] <- 1
        }
        maxf <- maxf + 1
        # print(maxf)
      }
      i <- 0
      while (i < N) {
        psx <- rt(1,1,mean[1]) # create a potential sample for a point in x
        psy <- rt(1,1,mean[2]) # create a potential sample for a point in y
        testsamp <- runif(1,0,suppressWarnings(maxf*dt(psx,1,mean[1])*dt(psy,1,mean[2]))) # create a test sample
        if (testsamp < f(c(psx,psy))) { # check that the test sample is below the pdf at the potential sample
          samples[i+1,1] <- psx
          samples[i+1,2] <- psy
          i <- i + 1
        }
      }
    }
    samples # outputs the samples matrix.
  }
}
