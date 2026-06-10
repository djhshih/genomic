suppressPackageStartupMessages(library(DNAcopy))

write_vec <- function(path, x) {
  write.table(data.frame(i=seq_along(x), value=x), file=path, sep="\t", row.names=FALSE, quote=FALSE)
}

run_case <- function(x, wts=NULL, alpha=0.01, nperm=200, min.width=2, kmax=25, nmin=200, eta=0.05, sbdry=FALSE) {
  chrom <- rep(1L, length(x))
  maploc <- seq_along(x)
  cna <- CNA(genomdat=x, chrom=chrom, maploc=maploc, data.type="logratio")
  if (is.null(wts)) {
    seg <- segment(cna, alpha=alpha, nperm=nperm, p.method="perm", min.width=min.width, kmax=kmax, nmin=nmin, eta=eta, sbdry=sbdry, verbose=0)
  } else {
    seg <- segment(cna, weights=wts, alpha=alpha, nperm=nperm, p.method="perm", min.width=min.width, kmax=kmax, nmin=nmin, eta=eta, sbdry=sbdry, verbose=0)
  }
  out <- seg$output
  data.frame(start=out$loc.start, end=out$loc.end, mean=out$seg.mean)
}

set.seed(1)
x3 <- c(rnorm(30, 0, 0.15), rnorm(20, 1.1, 0.15), rnorm(25, -0.9, 0.15), rnorm(25, 0.2, 0.15))
set.seed(2)
x4 <- c(rnorm(25, 0, 0.2), rnorm(25, 0.8, 0.2), rnorm(25, 0, 0.2), rnorm(25, 1.0, 0.2))

write_vec("tests/data/cbs_case3_noisy_input.tsv", x3)
write_vec("tests/data/cbs_case4_noisy_input.tsv", x4)
write.table(run_case(x3, NULL), file="tests/data/cbs_case3_noisy_expected.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(run_case(x4, NULL), file="tests/data/cbs_case4_noisy_expected.tsv", sep="\t", row.names=FALSE, quote=FALSE)
