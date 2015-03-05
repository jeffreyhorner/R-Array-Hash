## Get Jeffreys Image
FROM jeffreyhorner/r-array-hash

## Install requirements
RUN R -e 'install.packages(c("ggplot2", "microbenchmark", "jsonlite"))'
RUN RD -e 'install.packages(c("ggplot2", "microbenchmark", "jsonlite"))'

## The actual benchmarks
RUN R -e 'library(ggplot2); library(microbenchmark); library(jsonlite); print(microbenchmark(json <- toJSON(diamonds), fromJSON(json), times = 20));'
RUN RD -e 'library(ggplot2); library(microbenchmark); library(jsonlite); print(microbenchmark(json <- toJSON(diamonds), fromJSON(json), times = 20));'

