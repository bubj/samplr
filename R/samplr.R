samplr <- function(fun, a, b, N, c = NULL, d = NULL) {
  if (is.null(c) & is.null(d)) {
    g <- parse(text = fun) # changes the input fun string into an expression
    f <- function(x) {
      ifelse(a <= x & x <= b, eval(g[[1]]), 0) # defines the pdf as a function by evaluating the expression g
    }
    maxf <- optimize(f,c(a,b), tol = 0.0001, maximum = TRUE) # uses the optimize function to find the maximum of the pdf
    samples <- rep(0,N) # creates a vector for storing the sample values
    i <- 0
    while ( i < N) { # while loop that runs until N samples are selected
      potsamp <- runif(1,a,b) # creating a potential sample
      testsamp <- runif(1, 0, maxf$objective + .1) # creating a test sample. Uses maxf$objective to get the max of the function, and adds .1 to account for errors in tolerance
      if ( testsamp < f(potsamp) ) { # tests if the potential sample should be rejected. If it is kept, it is kept in the samples vector
        samples[i+1] = potsamp
        i = i + 1
      }
    }
    samples # output the samples
  } else {
    g <- parse(text = fun) # convert the input fun string into an expression
    f <- function(z) {
      x <- z[1]
      y <- z[2]
      ifelse(a <= x & x <= b & c <= y & y <= d, eval(g[[1]]), 0) # defining the pdf as a function where the input is a vector. The expression g is evaluated to get the equation.
    }
    maxf <- optim(c(b/2,d/2),f,control = list(fnscale = -1)) # finds the maximum of the pdf
    samples <- matrix(rep(0,2*N), nrow = N, ncol = 2) # creating a matrix to store the samples in
    i <- 0
    while( i < N) { # while loop that runs unil N accepted samples are found.
      psx <- runif(1,a,b) # creating a potential x value to be tested
      psy <- runif(1,c,d) # creating a potential y value to be tested
      testsamp <- runif(1,0, maxf$value + .1) # creating a test sample using maxf$value to find the max of the pdf and adding .1 to account for errors in tolerance.
      if (testsamp < f(c(psx,psy))) { # tests if the potential x and y values should be rejected. If kept, the x value is stored in the first column and the y value in the second of the samples matrix.
        samples[i+1,1] <- psx
        samples[i+1,2] <- psy
        i <- i +1
      }
    }
    samples # outputs the samples matrix.
  }
}
