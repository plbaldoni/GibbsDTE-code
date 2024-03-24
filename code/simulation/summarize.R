OK <- requireNamespace(c('devtools','BiocParallel'), quietly = TRUE)
if (!OK) {
  stop("devtools and BiocParallel packages required but is not installed (or can't be loaded)")
}
library(BiocParallel)
library(devtools)

load_all("../pkg")

workers <- 16
BPPARAM <- MulticoreParam(workers = workers,progressbar = TRUE)
register(BPPARAM)

summarizeSimulation(path = '../../output/simulation/data/',
                    dest = '../../output/simulation/summary/',
                    BPPARAM = BPPARAM)
