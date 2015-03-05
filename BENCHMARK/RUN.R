#
#
progname <- commandArgs()[1] # full path to R.

args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=4) stop("Need args: dataset hashsize run resultfile")

datafile <- args[1]
hashsize <- as.integer(args[2])
run      <- as.integer(args[3])
resfile  <- args[4]                  # will create if not exists

if (!file.exists(datafile)) stop("No file!")

# For calculating memory cache misses
#library(pbdPAPI)

set.seed(1)

# Test value for each hash entry
val <- paste(sample(letters,size=7),collapse='')


# The environment as hash to test
e <- new.env(hash=TRUE, size=hashsize)

constructHash <- function(filename, hashsize=2^10){
  con <- file(filename,open="r")
  on.exit(close(con))
  sum <- 0.0
  while(TRUE){
    line <- readLines(con,n=1)
    if (length(line) < 1) 
      break

    t1 <- proc.time()['elapsed']
    assign(line,val,envir=e);
    t2 <- proc.time()['elapsed']
    sum <- sum + (t2 - t1)
  }
  
  invisible(sum)
}

miss <- 0

searchHash <- function(filename){
  con <- file(filename,open="r")
  on.exit({close(con); lastline <<- line;})
  sum <- 0.0
  while(TRUE){
    line <- readLines(con,n=1)
    if (length(line) < 1) 
      break

    t1 <- proc.time()['elapsed']
    m <- get(line,envir=e)
    t2 <- proc.time()['elapsed']
    sum <- sum + (t2 - t1)

    # Miss should never miss!
    if (m!=val) miss <<- miss + 1
  }
 
  invisible(sum);
}

# Use this to see how many gc occurred
# gcinfo(TRUE)

invisible(gc(reset=TRUE)) # Reset is important here to estimate
                          # the hash memory size on next gc().

t1 <- proc.time()['elapsed']
x <- constructHash(datafile,hashsize)
t2 <- proc.time()['elapsed']
runsum <- t2 - t1;

z <- gc() # Estimation of hash memory size is z[11] + z[12]

t1 <- proc.time()['elapsed']
y <- searchHash(datafile)
t2 <- proc.time()['elapsed']
runsum <- runsum + (t2 - t1)

if (miss > 0) warning("WARNING!!!",miss,"missed variable lookups occurred!")

cachemiss <- 0
# Cache misses
#w <- system.cache(try(get('the',envir=e),silent=TRUE),gcFirst=TRUE,burnin=FALSE)
#cachemiss <- sum(unlist(w))

results <- data.frame(
  progname=progname,
  constructtime=x,
  searchtime=y,
  memory=z[11] + z[12],
  runtime=runsum,
  hashsize=hashsize,
  cachemiss=cachemiss,
  datafile=datafile,
  run=run
)

if (!file.exists(resfile)) {
  write.table(results,file=resfile,row.names=FALSE)
} else {
  write.table(results,file=resfile,col.names=FALSE,row.names=FALSE,append=TRUE)
}
