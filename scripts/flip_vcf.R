library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
in_vcf <- args[1] # grep -v '^##'
out_vcf <- args[2] # Need to concat back with all meta lines and col headers

vcf = read.table(in_vcf, header=TRUE, comment = "")

temp_col <- vcf$REF
vcf$REF <- vcf$ALT
vcf$ALT <- temp_col

swap_genotypes <- function(genotype) {
  case_when(
    genotype == '0/0' ~ '1/1',
    genotype == '1/1' ~ '0/0',
    TRUE ~ genotype
  )
}

num_columns <- ncol(vcf)
genotype_columns <- 10:(num_columns)
vcf[, genotype_columns] <- lapply(vcf[, genotype_columns], swap_genotypes)

write.table(vcf, file = out_vcf, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")




