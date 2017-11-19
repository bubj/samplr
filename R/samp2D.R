samp2D <- function(fun,a,b,c,d,N) {
  g <- parse(text = fun)
  f <- function(z) {
    x <- z[1]
    y <- z[2]
    ifelse(a <= x & x <= b & c <= y & y <= d, eval(g[[1]]), 0)
  }
  maxf <- optim(c(b/2,d/2),f,control = list(fnscale = -1))
  samples <- matrix(rep(0,2*N), nrow = N, ncol = 2)
  i <- 0
  while( i < N) {
    psx <- runif(1,a,b)
    psy <- runif(1,c,d)
    testsamp <- runif(1,0, maxf$value + .1)
    if (testsamp < f(c(psx,psy))) {
      samples[i+1,1] <- psx
      samples[i+1,2] <- psy
      i <- i +1
    }
  }
  samples
}
