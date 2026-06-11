suppressPackageStartupMessages(library(DNAcopy))

# Unified CBS fixture generator.
# Generates deterministic input/expected TSVs for multiple parameterizations
# used by tests/cbs_test.cpp.

write_vec <- function(path, x, value.name="value") {
  write.table(data.frame(index=seq_along(x), setNames(list(x), value.name)),
              file=path, sep="\t", row.names=FALSE, quote=FALSE)
}

run_case <- function(x,
                     wts=NULL,
                     alpha=0.01,
                     nperm=200,
                     p.method="perm",
                     min.width=2,
                     kmax=25,
                     nmin=200,
                     eta=0.05,
                     sbdry=FALSE) {
  chrom <- rep(1L, length(x))
  maploc <- seq_along(x)
  cna <- CNA(genomdat=x, chrom=chrom, maploc=maploc, data.type="logratio")
  seg <- if (is.null(wts)) {
    segment(cna,
            alpha=alpha,
            nperm=nperm,
            p.method=p.method,
            min.width=min.width,
            kmax=kmax,
            nmin=nmin,
            eta=eta,
            sbdry=sbdry,
            verbose=0)
  } else {
    segment(cna,
            weights=wts,
            alpha=alpha,
            nperm=nperm,
            p.method=p.method,
            min.width=min.width,
            kmax=kmax,
            nmin=nmin,
            eta=eta,
            sbdry=sbdry,
            verbose=0)
  }
  out <- seg$output
  data.frame(start=out$loc.start, end=out$loc.end, mean=out$seg.mean)
}

write_expected <- function(path, df) {
  write.table(df, file=path, sep="\t", row.names=FALSE, quote=FALSE)
}

write_case_bundle <- function(name, x, wts=NULL, configs) {
  write_vec(file.path("tests/data", paste0(name, "_input.tsv")), x)
  if (!is.null(wts)) {
    write_vec(file.path("tests/data", paste0(name, "_weights.tsv")), wts, value.name="weight")
  }
  for (cfg in configs) {
    suffix <- cfg$suffix
    out <- run_case(x,
                    wts=wts,
                    alpha=cfg$alpha,
                    nperm=cfg$nperm,
                    p.method=cfg$p.method,
                    min.width=cfg$min.width,
                    kmax=cfg$kmax,
                    nmin=cfg$nmin,
                    eta=cfg$eta,
                    sbdry=cfg$sbdry)
    write_expected(file.path("tests/data", paste0(name, suffix, "_expected.tsv")), out)
    message("Generated ", name, suffix, "_expected.tsv")
    print(out)
  }
}

perm_cfg <- list(suffix="", alpha=0.01, nperm=200, p.method="perm", min.width=2, kmax=25, nmin=200, eta=0.05, sbdry=FALSE)
perm_alt_cfg <- list(suffix="_perm_alt", alpha=0.05, nperm=100, p.method="perm", min.width=3, kmax=25, nmin=200, eta=0.05, sbdry=FALSE)
hybrid_cfg <- list(suffix="_hybrid", alpha=0.01, nperm=200, p.method="hybrid", min.width=2, kmax=25, nmin=200, eta=0.05, sbdry=FALSE)
hybrid_alt_cfg <- list(suffix="_hybrid_alt", alpha=0.05, nperm=100, p.method="hybrid", min.width=3, kmax=25, nmin=200, eta=0.05, sbdry=FALSE)
noisy_sbdry <- rep(201L, 201 * 202 / 2 + 2)
noisy_perm_cfg <- list(suffix="", alpha=0.01, nperm=200, p.method="perm", min.width=2, kmax=25, nmin=200, eta=0.05, sbdry=noisy_sbdry)
noisy_perm_alt_cfg <- list(suffix="_perm_alt", alpha=0.05, nperm=100, p.method="perm", min.width=3, kmax=25, nmin=200, eta=0.05, sbdry=noisy_sbdry)
noisy_hybrid_cfg <- list(suffix="_hybrid", alpha=0.01, nperm=200, p.method="hybrid", min.width=2, kmax=25, nmin=200, eta=0.05, sbdry=noisy_sbdry)
noisy_hybrid_alt_cfg <- list(suffix="_hybrid_alt", alpha=0.05, nperm=100, p.method="hybrid", min.width=3, kmax=25, nmin=200, eta=0.05, sbdry=noisy_sbdry)

x1 <- c(rep(0.0, 20), rep(1.5, 20), rep(0.0, 20))
x2 <- c(rep(0.0, 15), rep(2.0, 15), rep(-1.5, 15), rep(0.0, 15))
w2 <- c(rep(1.0, 15), rep(0.5, 15), rep(2.0, 15), rep(1.0, 15))

set.seed(1)
x3 <- c(rnorm(30, 0, 0.15), rnorm(20, 1.1, 0.15), rnorm(25, -0.9, 0.15), rnorm(25, 0.2, 0.15))
set.seed(2)
x4 <- c(rnorm(25, 0, 0.2), rnorm(25, 0.8, 0.2), rnorm(25, 0, 0.2), rnorm(25, 1.0, 0.2))

write_case_bundle("cbs_case1", x1, NULL, list(perm_cfg, perm_alt_cfg, hybrid_cfg, hybrid_alt_cfg))
write_case_bundle("cbs_case2_weighted", x2, w2, list(perm_cfg, perm_alt_cfg, hybrid_cfg, hybrid_alt_cfg))
write_case_bundle("cbs_case3_noisy", x3, NULL, list(noisy_perm_cfg, noisy_perm_alt_cfg, noisy_hybrid_cfg, noisy_hybrid_alt_cfg))
write_case_bundle("cbs_case4_noisy", x4, NULL, list(noisy_perm_cfg, noisy_perm_alt_cfg, noisy_hybrid_cfg, noisy_hybrid_alt_cfg))
