suppressPackageStartupMessages(library(DNAcopy))

write_raw_matrix <- function(path, marker, chrom_name, chrom_num, position, samples) {
  df <- data.frame(marker=marker, chromosome=chrom_name, position=position, check.names=FALSE)
  for (nm in names(samples)) df[[nm]] <- samples[[nm]]
  write.table(df, file=path, sep="\t", row.names=FALSE, quote=FALSE)
}

write_seg_output <- function(path, rows) {
  write.table(rows, file=path, sep="\t", row.names=FALSE, quote=FALSE)
}

run_sample <- function(values, chrom, maploc,
                       smooth.region=10, outlier.SD.scale=4, smooth.SD.scale=2, trim=0.025,
                       alpha=0.01, nperm=200, p.method="perm", min.width=2, kmax=25, nmin=200, eta=0.05) {
  cna <- CNA(genomdat=values, chrom=chrom, maploc=maploc, data.type="logratio", presorted=TRUE)
  sm <- smooth.CNA(cna,
                   smooth.region=smooth.region,
                   outlier.SD.scale=outlier.SD.scale,
                   smooth.SD.scale=smooth.SD.scale,
                   trim=trim)
  seg <- segment(sm,
                 alpha=alpha,
                 nperm=nperm,
                 p.method=p.method,
                 min.width=min.width,
                 kmax=kmax,
                 nmin=nmin,
                 eta=eta,
                 verbose=0)
  out <- seg$output
  data.frame(
    ID=out$ID,
    chrom=out$chrom,
    loc.start=out$loc.start,
    loc.end=out$loc.end,
    num.mark=out$num.mark,
    seg.mean=out$seg.mean,
    check.names=FALSE
  )
}

chrom <- c(rep(1L, 8), rep(2L, 8))
chrom_name <- c(rep("chr1", 8), rep("chr2", 8))
position <- c(seq(10, 80, by=10), seq(10, 80, by=10))
marker <- sprintf("m%02d", seq_along(chrom))

sample1 <- c(-0.2, -0.1, 0.0, 1.1, 1.2, 1.0, 0.1, -0.1,
             -0.3, -0.2, 0.0, 0.9, 1.0, 1.1, 0.2, -0.2)
sample2 <- c(-0.4, -0.3, -0.2, 0.0, 0.1, 0.2, 1.4, 1.5,
             -0.5, -0.4, -0.2, 0.0, 0.2, 0.1, 1.3, 1.4)

input_path <- "tests/data/segment_cli_case1_input.cn"
expected_path <- "tests/data/segment_cli_case1_expected.seg"

write_raw_matrix(input_path, marker, chrom_name, chrom, position, list(sample1=sample1, sample2=sample2))

expected <- rbind(
  run_sample(sample1, chrom, position),
  run_sample(sample2, chrom, position)
)

write_seg_output(expected_path, expected)
