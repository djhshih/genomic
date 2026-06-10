suppressPackageStartupMessages(library(DNAcopy))

run_case <- function(name, x, wts=NULL, alpha=0.01, nperm=200, min.width=2, kmax=25, nmin=200, eta=0.05, sbdry=FALSE) {
  chrom <- rep(1L, length(x))
  maploc <- seq_along(x)
  cna <- CNA(genomdat=x, chrom=chrom, maploc=maploc, data.type="logratio")
  if (is.null(wts)) {
    seg <- segment(cna,
                   alpha=alpha,
                   nperm=nperm,
                   p.method="perm",
                   min.width=min.width,
                   kmax=kmax,
                   nmin=nmin,
                   eta=eta,
                   sbdry=sbdry,
                   verbose=0)
  } else {
    seg <- segment(cna,
                   weights=wts,
                   alpha=alpha,
                   nperm=nperm,
                   p.method="perm",
                   min.width=min.width,
                   kmax=kmax,
                   nmin=nmin,
                   eta=eta,
                   sbdry=sbdry,
                   verbose=0)
  }
  out <- seg$output
  write.table(data.frame(index=seq_along(x), value=x),
              file=paste0("tests/data/", name, "_input.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  if (!is.null(wts)) {
    write.table(data.frame(index=seq_along(wts), weight=wts),
                file=paste0("tests/data/", name, "_weights.tsv"),
                sep="\t", row.names=FALSE, quote=FALSE)
  }
  write.table(data.frame(start=out$loc.start, end=out$loc.end, mean=out$seg.mean),
              file=paste0("tests/data/", name, "_expected.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
}

set.seed(1)
x1 <- c(rep(0.0, 20), rep(1.5, 20), rep(0.0, 20))
x2 <- c(rep(0.0, 15), rep(2.0, 15), rep(-1.5, 15), rep(0.0, 15))
w2 <- c(rep(1.0, 15), rep(0.5, 15), rep(2.0, 15), rep(1.0, 15))
x3 <- c(rep(0.0, 30), rep(1.2, 25), rep(-0.8, 20), rep(0.3, 25)) + rnorm(100, mean=0, sd=0.15)
x4 <- c(rep(-0.5, 20), rep(0.8, 20), rep(1.6, 20), rep(0.2, 20), rep(-1.0, 20)) + rnorm(100, mean=0, sd=0.2)

run_case("cbs_case1", x1, NULL)
run_case("cbs_case2_weighted", x2, w2)
run_case("cbs_case3_noisy", x3, NULL)
run_case("cbs_case4_noisy", x4, NULL)
