OK <- requireNamespace("devtools", quietly = TRUE)
if (!OK) {
  stop("devtools package required but is not installed (or can't be loaded)")
}

################################################################################
# Getting parameters
args <- R.utils::commandArgs(asValues = TRUE)
print(args)
genome <- args[['genome']]
read.length <- as.integer(args[['rlen']])
fc <- as.numeric(args[['fc']])
paired.end <- as.logical(args[['pe']])
max.tx <- as.integer(args[['mtx']])
scenario <- args[['scenario']]
libs.per.group <- as.integer(args[['libs']])
workers <- as.integer(args[['workers']])
projdir <- normalizePath(file.path(dirname(as.character(args[['file']])),"../.."))
################################################################################

tmpdir <- tempdir(check = TRUE)
print(tmpdir)

devtools::load_all(file.path(projdir,"code/pkg"))

if (genome == 'mm39') {
  fasta <- file.path(projdir,'data/annotation/mm39/gencode.vM27.transcripts.fa.gz')
} else {
  stop('only simulations for mm39 were run')
}

n.libs <- rep(libs.per.group,2)

if (scenario == 'balanced') {
  lib.sizes <- rep(50e6,sum(n.libs))
}
if (scenario == 'unbalanced') {
  lib.sizes <- rep(rep(c(25e6,100e6),length.out = libs.per.group),2)
}

dest <- normalizePath('.')
dir.create(dest,showWarnings = FALSE,recursive = TRUE)

seed <- sample(-999999999:999999999,1)
message("Seed: ",seed)

simulateExperiment(dest = dest,
                   fasta = fasta,
                   genome = genome,
                   tmpdir = tmpdir,
                   max.tx = max.tx,
                   workers = workers,
                   fc = fc,
                   n.libs = n.libs,
                   lib.sizes = lib.sizes,
                   paired.end = paired.end,
                   keep.fastq = FALSE,
                   read.length = read.length,
                   top.cutoff = 30000,
                   num.DE = 3000,
                   seed = seed,
                   run.kallisto = FALSE,
                   run.salmon = TRUE,
                   bcv.true = 0.6,
                   seq.b = c(2,3),
                   opts.salmon.boot =  paste('-p',workers,'-l A --numBootstraps 20 --validateMappings'),
                   opts.salmon.gibbs =  paste('-p',workers,'-l A --numGibbsSamples 20 --validateMappings'),
                   bin.salmon = "/stornext/General/data/academic/lab_smyth/baldoni.p/software/salmon-1.10.0_linux_x86_64/bin/salmon",
                   index.salmon = "/stornext/General/data/academic/lab_smyth/baldoni.p/software/SalmonIndex/mm39-M27/transcripts_index",
                   bin.kallisto = "/stornext/General/data/academic/lab_smyth/baldoni.p/software/kallisto/kallisto",
                   index.kallisto = "/stornext/General/data/academic/lab_smyth/baldoni.p/software/kallistoIndex/mm39/transcripts_index")
