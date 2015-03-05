## Get Jeffreys Image
FROM jeffreyhorner/r-array-hash

## Install requirements
RUN mkdir -p ~/R312 && \
    mkdir -p ~/R320 && \
    R -e '.libPaths("~/R312"); install.packages(c("ggplot2", "microbenchmark", "jsonlite"))' && \
    RD -e '.libPaths("~/R320"); install.packages(c("ggplot2", "microbenchmark", "jsonlite"))'

## The actual benchmarks
RUN R -e '.libPaths("~/R312"); library(ggplot2); library(microbenchmark); library(jsonlite); json <- toJSON(diamonds); print(microbenchmark(toJSON(diamonds), fromJSON(json), times = 20));' && \
    RD -e '.libPaths("~/R320"); library(ggplot2); library(microbenchmark); library(jsonlite); json <- toJSON(diamonds); print(microbenchmark(toJSON(diamonds), fromJSON(json), times = 20));'
