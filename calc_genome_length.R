# Calculate genome length from a .bim file
# Author: SSG
#------------
# This script accepts a PLINK .bim file and returns the number of bases that the file spans
# The required argument is: Input (name of file)
# Ex.Rscript calc_genome_length.R dat_autosomes_HWE-MAFfiltered.bim

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Update defaults, check for correct number of arguments
if (length(args)<1) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (length(args)>1) {
  stop("Too many arguments", call.=FALSE)
}

options(stringsAsFactors=F)
bim <- read.table(args[1])

genome.length <- 0
for (chr in unique(bim$V1)) {
  chr.bim <- bim[bim$V1==chr,]
    chr.bim <- chr.bim[order(chr.bim$V4),]
    genome.length <- genome.length + (chr.bim[nrow(chr.bim),4] - chr.bim[1,4])
}

cat(genome.length, "\n")
