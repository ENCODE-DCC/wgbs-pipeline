library(bsseq)
library(optparse)

options = list(
  make_option(c("-i", "--infile"),
              help = "input BED file name in Bismark coverage format",),
  make_option(c("-o", "--outfile"),
              help = "name of file to output smoothed methylation estimate to",),
  make_option(
    c("-w", "--numworkers"),
    type = "integer",
    default = 8,
    help = "Number of BiocParallelBackend workers to use [default %default]",
  ),
  make_option(
    c("-t", "--numthreads"),
    type = "integer",
    default = 2,
    help = "Number of threads per worker to use [default %default]",
  )
)

arg_parser <- OptionParser(option_list = options)
args <- parse_args(arg_parser)

bsseq <-
  read.bismark(
    files = args$infile,
    colData = DataFrame(row.names = "test_data"),
    rmZeroCov = FALSE,
    strandCollapse = FALSE,
    verbose = TRUE,
    nThread = args$numthreads,
    BPPARAM = MulticoreParam(workers = args$numworkers, progressbar = TRUE),
    BACKEND = "HDF5Array"
  )

# Seems to use either the minimum of NUM_CORES cores or the number of available CPUs on
# the machine. Fits the BSmooth model
bsseq.fit <- BSmooth(bsseq,
                     BPPARAM = MulticoreParam(workers = args$numworkers, progressbar = TRUE),)

# Obtain the smoothed methylation profile. This will replace any positions with NaNs
# (where there was no coverage in input bed). Conceptually similar to sklearn transform
# methods
meth <- getMeth(bsseq.fit, type = "smooth", what = "perBase")

# meth is a DelayedMatrix, need to materialize in memory to write to a TSV. See
# documentation for DelayedArray-class,  "In-memory versus on-disk realization", p. 13
# https://bioconductor.org/packages/devel/bioc/manuals/DelayedArray/man/DelayedArray.pdf
write.table(
  as.array(meth),
  file = args$outfile,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
)
