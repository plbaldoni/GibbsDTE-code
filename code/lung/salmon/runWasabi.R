OK <- requireNamespace("wasabi", quietly = TRUE)
if (!OK) {
  stop("wasabi package required but is not installed (or can't be loaded). Download the package from https://github.com/COMBINE-lab/wasabi")
}
library(wasabi)
prepare_fish_for_sleuth(list.dirs('../../../output/lung/salmon',recursive = FALSE))
sessionInfo()
