library(Rsubread)
buildindex(basename = "../../../data/lung/index-subread/GRCh38.genome_index",
           reference = "../../../data/lung/index/GRCh38.p13.genome_sequins.fa.gz",memory = 64000)
sessionInfo()
