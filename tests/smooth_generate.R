suppressPackageStartupMessages(library(DNAcopy))

write_vec <- function(path, x, value.name="value") {
  write.table(data.frame(index=seq_along(x), setNames(list(x), value.name)),
              file=path, sep="\t", row.names=FALSE, quote=FALSE)
}

write_matrix <- function(path, cols) {
  df <- data.frame(index=seq_along(cols[[1]]))
  for (nm in names(cols)) df[[nm]] <- cols[[nm]]
  write.table(df, file=path, sep="\t", row.names=FALSE, quote=FALSE)
}

run_smooth <- function(values, chrom, smooth.region=10, outlier.SD.scale=4,
                       smooth.SD.scale=2, trim=0.025) {
  maploc <- seq_along(values)
  cna <- CNA(genomdat=values, chrom=chrom, maploc=maploc, data.type="logratio", presorted=TRUE)
  as.numeric(smooth.CNA(cna,
                        smooth.region=smooth.region,
                        outlier.SD.scale=outlier.SD.scale,
                        smooth.SD.scale=smooth.SD.scale,
                        trim=trim)[,3])
}

run_smooth_matrix <- function(samples, chrom, smooth.region=10, outlier.SD.scale=4,
                              smooth.SD.scale=2, trim=0.025) {
  maploc <- seq_along(chrom)
  cna <- CNA(genomdat=as.matrix(samples), chrom=chrom, maploc=maploc, data.type="logratio", presorted=TRUE)
  out <- smooth.CNA(cna,
                    smooth.region=smooth.region,
                    outlier.SD.scale=outlier.SD.scale,
                    smooth.SD.scale=smooth.SD.scale,
                    trim=trim)
  as.data.frame(out[, -(1:2), drop=FALSE])
}

chrom1 <- c(rep(1L, 8), rep(2L, 6))
x1 <- c(0, 0.05, -0.04, 5, 0.02, -0.03, 0.01, 0.00, 0, 0.02, 0.01, -4.5, -0.02, 0.00)
write_vec("tests/data/smooth_case1_input.tsv", x1)
write_vec("tests/data/smooth_case1_chrom.tsv", chrom1, "chrom")
write_vec("tests/data/smooth_case1_expected.tsv", run_smooth(x1, chrom1), "value")

chrom2 <- c(rep(1L, 5), rep(2L, 5))
x2 <- c(0, NA, 0.1, 6, 0.0, 0.2, Inf, 0.1, -5, 0.0)
write_vec("tests/data/smooth_case2_input.tsv", x2)
write_vec("tests/data/smooth_case2_chrom.tsv", chrom2, "chrom")
write_vec("tests/data/smooth_case2_expected.tsv", run_smooth(x2, chrom2, smooth.region=2), "value")

chrom3 <- c(rep(1L, 7), rep(2L, 7))
s1 <- c(0, 0, 0.1, 4.5, 0.0, 0.1, 0.0, 0, 0, -3.5, 0.1, 0.0, 0.2, 0.0)
s2 <- c(0.2, 0.0, -0.1, 5.0, 0.1, -0.1, 0.0, 0.1, 0.0, -4.0, 0.0, 0.2, 0.1, 0.0)
write_matrix("tests/data/smooth_case3_input.tsv", list(sample1=s1, sample2=s2))
write_vec("tests/data/smooth_case3_chrom.tsv", chrom3, "chrom")
expected3 <- run_smooth_matrix(data.frame(sample1=s1, sample2=s2), chrom3, smooth.region=2)
write_matrix("tests/data/smooth_case3_expected.tsv", list(sample1=expected3[[1]], sample2=expected3[[2]]))
