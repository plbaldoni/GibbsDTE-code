OK <- requireNamespace("data.table", quietly = TRUE)
if (!OK) {
  stop("data.table package required but is not installed (or can't be loaded)")
}
library(data.table)

dest.data <- '../../output/simulation-large/data'

for (genome in c('mm39')) {
  for (len in c(100)) {
    for (fc in c(2)) {
      for (scenario in c('unbalanced')) {
        for (paired.end in c(TRUE)) {
          for (max.tx in c(9999)) {
            for (libs.per.group in c(100)) {
              for (simulation in 1:20) {
                dest <- file.path(dest.data,
                                  genome,
                                  paste0('readlen-',len),
                                  paste0('fc',gsub("\\.","_",as.character(fc))),
                                  paste0(ifelse(paired.end,'paired','single'),'-end'),
                                  paste0(max.tx,'TxPerGene'),
                                  scenario,
                                  paste0(libs.per.group,'libsPerGroup'),
                                  paste0('simulation-',simulation))
                dir.create(dest,recursive = TRUE,showWarnings = FALSE)

                if (!length(list.files(dest, 'time.tsv', recursive = TRUE)) == 6) {
                  dt <- data.table(dest = dest,
                                   genome = genome,
                                   rlen = len,
                                   fc = fc,
                                   pe = paired.end,
                                   mtx = max.tx,
                                   scenario = scenario,
                                   libs = libs.per.group,
                                   simulation = simulation)

                  fwrite(x = dt,
                         file = 'parameters.txt',
                         sep = "\t",
                         quote = FALSE,
                         append = TRUE,
                         row.names = FALSE)
                }
              }
            }
          }
        }
      }
    }
  }
}
