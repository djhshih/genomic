summarize_cn_at_position_local <- function(seg, pos, direction, cutoff) {
  idx <- seg$start <= pos & pos <= seg$end
  logr <- direction * seg$logr[idx]
  idx2 <- logr > cutoff
  if (sum(idx2) > 0) {
    sum(exp(logr[idx2])) / length(logr)
  } else {
    0
  }
}

summarize_cn_local <- function(seg, direction, cutoff, positions=NULL) {
  if (is.null(positions)) {
    positions <- sort(unique(c(seg$start, seg$end)))
  }
  values <- unlist(lapply(positions, function(pos) {
    summarize_cn_at_position_local(seg, pos, direction, cutoff)
  }))
  data.frame(position=positions, value=values)
}

write_seg_like <- function(path, df) {
  write.table(df, file=path, sep="\t", row.names=FALSE, quote=FALSE)
}

write_summary <- function(path, d) {
  write.table(d, file=path, sep="\t", row.names=FALSE, quote=FALSE)
}

seg_case1 <- data.frame(
  sample=c("s1","s1","s1"),
  chromosome=c("1","1","1"),
  start=c(10,20,30),
  end=c(25,35,40),
  nprobes=c(5,5,4),
  logr=c(0.2,0.8,0.1)
)
write_seg_like("tests/data/cngpld_case1_input.seg", seg_case1)
write_summary("tests/data/cngpld_case1_amp_expected.tsv", summarize_cn_local(seg_case1, direction=1, cutoff=0.5))
write_summary("tests/data/cngpld_case1_del_expected.tsv", summarize_cn_local(seg_case1, direction=-1, cutoff=0.5))

seg_case2 <- data.frame(
  sample=c("s1","s1","s1"),
  chromosome=c("1","1","1"),
  start=c(10,20,30),
  end=c(25,35,40),
  nprobes=c(5,5,4),
  logr=c(-0.2,-0.9,-0.1)
)
write_seg_like("tests/data/cngpld_case2_input.seg", seg_case2)
write_summary("tests/data/cngpld_case2_amp_expected.tsv", summarize_cn_local(seg_case2, direction=1, cutoff=0.5))
write_summary("tests/data/cngpld_case2_del_expected.tsv", summarize_cn_local(seg_case2, direction=-1, cutoff=0.5))

seg_case3 <- data.frame(
  sample=c("s1","s1","s1"),
  chromosome=c("1","1","1"),
  start=c(10,10,10),
  end=c(20,20,20),
  nprobes=c(3,3,3),
  logr=c(0.7,0.2,0.1)
)
write_seg_like("tests/data/cngpld_case3_input.seg", seg_case3)
write_summary("tests/data/cngpld_case3_amp_expected.tsv", summarize_cn_local(seg_case3, direction=1, cutoff=0.5, positions=c(10,20)))

seg_case4 <- data.frame(
  sample=c("s1","s1"),
  chromosome=c("1","1"),
  start=c(100,200),
  end=c(120,220),
  nprobes=c(2,2),
  logr=c(1.0,-1.0)
)
write_seg_like("tests/data/cngpld_case4_input.seg", seg_case4)
write_summary("tests/data/cngpld_case4_amp_expected.tsv", summarize_cn_local(seg_case4, direction=1, cutoff=0.5, positions=c(50,100,120,150,200,220,250)))
