gc(reset=TRUE)
#gcinfo(TRUE)
runs <- 3
times <- rep(0, 15); dim(times) <- c(5,3)
cumulate <- 0; b <- 0
for (i in 1:runs) {
  b <- numeric(1000);
  invisible(gc())
  timing <- system.time({
  	# Rem: there are faster ways to do this
  	# but here we want to time loops (220*220 'for' loops)! 
    for (j in 1:1000) {
      for (k in 1:1000) {
        b <- b + 1
      }
    }
  })[3]
  cumulate <- cumulate + timing
}
timing <- cumulate/runs
times[4, 3] <- timing
cat(c("Creation of a simple sum (loops)_______ (sec): ", timing, "\n"))
gc()
