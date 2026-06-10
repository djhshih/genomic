suppressPackageStartupMessages(library(IRanges))

write_case <- function(name, starts, ends, queries) {
  subject <- IRanges(start=starts, end=ends)
  input_df <- data.frame(id=seq_along(starts) - 1L, start=starts, end=ends)
  write.table(input_df,
              file=paste0("tests/data/", name, "_intervals.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)

  expected_qid <- integer()
  expected_hit <- integer()
  for (qi in seq_len(nrow(queries))) {
    query <- IRanges(start=queries$qstart[qi], end=queries$qend[qi])
    hits <- subjectHits(findOverlaps(query, subject, type="any", select="all")) - 1L
    if (length(hits) != 0L) {
      expected_qid <- c(expected_qid, rep(qi - 1L, length(hits)))
      expected_hit <- c(expected_hit, hits)
    }
  }
  expected <- data.frame(qid=expected_qid, hit=expected_hit)
  write.table(queries,
              file=paste0("tests/data/", name, "_queries.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  write.table(expected,
              file=paste0("tests/data/", name, "_expected.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
}

case1_queries <- data.frame(
  qstart=c(1L, 5L, 10L, 21L, 50L, 90L),
  qend=c(1L, 5L, 20L, 25L, 60L, 95L)
)
write_case(
  "nclist_case1",
  c(1L, 10L, 15L, 30L, 50L, 70L),
  c(100L, 20L, 25L, 30L, 60L, 90L),
  case1_queries
)

case2_queries <- data.frame(
  qstart=c(4L, 20L, 23L, 101L, 109L),
  qend=c(4L, 20L, 23L, 105L, 112L)
)
write_case(
  "nclist_case2",
  c(1L, 2L, 4L, 20L, 21L, 100L, 100L),
  c(10L, 3L, 5L, 30L, 22L, 110L, 100L),
  case2_queries
)

case3_queries <- data.frame(
  qstart=c(5L, 6L, 7L, 11L, 12L, 14L),
  qend=c(5L, 6L, 9L, 11L, 13L, 14L)
)
write_case(
  "nclist_case3",
  c(5L, 5L, 5L, 6L, 10L, 12L),
  c(5L, 10L, 15L, 6L, 12L, 14L),
  case3_queries
)
