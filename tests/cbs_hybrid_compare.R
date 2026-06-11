suppressPackageStartupMessages(library(DNAcopy))

run_case <- function(x,
                     wts=NULL,
                     alpha=0.01,
                     nperm=200,
                     min.width=2,
                     kmax=25,
                     nmin=200,
                     eta=0.05,
                     sbdry=FALSE) {
  chrom <- rep(1L, length(x))
  maploc <- seq_along(x)
  cna <- CNA(genomdat=x, chrom=chrom, maploc=maploc, data.type="logratio")
  seg <- if (is.null(wts)) {
    segment(cna, alpha=alpha, nperm=nperm, p.method="hybrid", min.width=min.width, kmax=kmax, nmin=nmin, eta=eta, sbdry=sbdry, verbose=0)
  } else {
    segment(cna, weights=wts, alpha=alpha, nperm=nperm, p.method="hybrid", min.width=min.width, kmax=kmax, nmin=nmin, eta=eta, sbdry=sbdry, verbose=0)
  }
  out <- seg$output
  data.frame(start=out$loc.start, end=out$loc.end, mean=out$seg.mean)
}

x1 <- read.table("tests/data/cbs_case1_input.tsv", header=TRUE, sep="\t")$value
x2 <- read.table("tests/data/cbs_case2_weighted_input.tsv", header=TRUE, sep="\t")$value
w2 <- read.table("tests/data/cbs_case2_weighted_weights.tsv", header=TRUE, sep="\t")$weight

cat("case1 hybrid\n")
print(run_case(x1, alpha=0.01, nperm=200, min.width=2, kmax=25, nmin=200, eta=0.05))
cat("case1 hybrid alt\n")
print(run_case(x1, alpha=0.05, nperm=100, min.width=3, kmax=25, nmin=200, eta=0.05))
cat("case2 weighted hybrid\n")
print(run_case(x2, w2, alpha=0.01, nperm=200, min.width=2, kmax=25, nmin=200, eta=0.05))
cat("case2 weighted hybrid alt\n")
print(run_case(x2, w2, alpha=0.05, nperm=100, min.width=3, kmax=25, nmin=200, eta=0.05))
