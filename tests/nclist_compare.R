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

case4_intervals <- data.frame(
  id=0:39,
  start=c(18L, 98L, 64L, 61L, 27L, 50L, 90L, 35L, 30L, 14L,
          4L, 70L, 49L, 55L, 68L, 98L, 64L, 30L, 4L, 84L,
          88L, 55L, 68L, 98L, 64L, 35L, 30L, 14L, 4L, 70L,
          49L, 55L, 68L, 98L, 64L, 30L, 4L, 84L, 88L, 12L),
  end=c(54L, 102L, 95L, 89L, 58L, 78L, 118L, 57L, 53L, 37L,
        20L, 105L, 83L, 75L, 100L, 120L, 97L, 46L, 11L, 101L,
        94L, 86L, 72L, 132L, 89L, 61L, 34L, 47L, 18L, 78L,
        60L, 84L, 92L, 125L, 88L, 66L, 22L, 110L, 115L, 41L)
)
case4_queries <- data.frame(
  qstart=c(1L, 7L, 12L, 18L, 24L, 30L, 36L, 42L, 48L, 54L,
           60L, 66L, 72L, 78L, 84L, 90L, 96L, 102L, 108L, 114L,
           120L, 5L, 33L, 57L, 81L, 99L, 117L, 15L, 45L, 75L),
  qend=c(6L, 13L, 19L, 27L, 31L, 39L, 45L, 50L, 59L, 63L,
         68L, 74L, 80L, 87L, 91L, 97L, 104L, 109L, 116L, 122L,
         128L, 25L, 44L, 70L, 95L, 118L, 130L, 29L, 56L, 90L)
)
write_case(
  "nclist_case4",
  case4_intervals$start,
  case4_intervals$end,
  case4_queries
)

case5_queries <- data.frame(
  qstart=c(1L, 2L, 3L, 4L, 5L, 6L, 10L, 20L, 40L, 80L, 120L),
  qend=c(1L, 2L, 3L, 4L, 5L, 6L, 10L, 20L, 40L, 80L, 120L)
)
write_case(
  "nclist_case5_deep",
  c(1L, 2L, 3L, 4L, 5L, 6L),
  c(120L, 110L, 100L, 90L, 80L, 70L),
  case5_queries
)

case6_queries <- data.frame(
  qstart=c(10L, 20L, 30L, 40L, 50L, 60L, 70L, 80L),
  qend=c(10L, 25L, 35L, 40L, 55L, 60L, 75L, 80L)
)
write_case(
  "nclist_case6_same_start",
  c(10L, 10L, 10L, 10L, 10L, 10L),
  c(15L, 25L, 35L, 45L, 60L, 80L),
  case6_queries
)

case7_queries <- data.frame(
  qstart=c(1L, 10L, 20L, 30L, 40L, 50L, 60L, 70L),
  qend=c(10L, 20L, 30L, 40L, 50L, 60L, 70L, 80L)
)
write_case(
  "nclist_case7_same_end",
  c(1L, 11L, 21L, 31L, 41L, 51L, 61L, 71L),
  c(80L, 80L, 80L, 80L, 80L, 80L, 80L, 80L),
  case7_queries
)

case8_queries <- data.frame(
  qstart=c(5L, 15L, 25L, 35L, 45L, 55L, 65L, 75L, 85L),
  qend=c(8L, 18L, 28L, 38L, 48L, 58L, 68L, 78L, 88L)
)
write_case(
  "nclist_case8_sparse",
  c(1L, 20L, 40L, 60L, 80L),
  c(10L, 30L, 50L, 70L, 90L),
  case8_queries
)

case9_queries <- data.frame(
  qstart=c(1L, 8L, 12L, 18L, 24L, 32L, 38L),
  qend=c(5L, 14L, 20L, 26L, 30L, 36L, 42L)
)
write_case(
  "nclist_case9_dense",
  c(1L, 5L, 9L, 13L, 17L, 21L, 25L, 29L, 33L),
  c(10L, 14L, 18L, 22L, 26L, 30L, 34L, 38L, 42L),
  case9_queries
)

set.seed(42)
for (case_idx in 1:3) {
  n_intervals <- 120L
  n_queries <- 80L
  starts <- sample(1:400, n_intervals, replace=TRUE)
  widths <- sample(0:40, n_intervals, replace=TRUE)
  ends <- starts + widths
  qstarts <- sample(1:400, n_queries, replace=TRUE)
  qwidths <- sample(0:30, n_queries, replace=TRUE)
  qends <- qstarts + qwidths
  write_case(
    paste0("nclist_fuzz", case_idx),
    starts,
    ends,
    data.frame(qstart=qstarts, qend=qends)
  )
}
