x <- read.table("genome-info_snp6_na24_hg18_combined_cnv-filtered.txt", sep="\t", as.is=TRUE);

x1 <- x[seq(1, nrow(x), 10000), ];
f1 <- x[sample(1:nrow(x1), 50), ];
y1 <- x1[!(x1$marker %in% f1$marker),];
y1$position <- as.numeric(y1$position);
y1$chromosome <- factor(y1$chromosome, levels=c(1:22, "X", "Y", "MT"));
y1 <- y1[order(y1$chromosome, y1$position), ];
write.table(x1, "cnfilt1a.in", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);
write.table(f1, "cnfilt1b.in", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);
write.table(y1, "cnfilt1.ans", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);


x2 <- x[sample(1:nrow(x), 200),]
f2 <- x2[sample(1:nrow(x2), 100),];
y2 <- x2[!(x2$marker %in% f2$marker),];
y2$position <- as.numeric(y2$position);
y2$chromosome <- factor(y2$chromosome, levels=c(1:22, "X", "Y", "MT"));
y2 <- y2[order(y2$chromosome, y2$position), ];
write.table(x2, "cnfilt2a.in", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);
write.table(f2, "cnfilt2b.in", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);
write.table(y2, "cnfilt2.ans", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);

x3 <- x[sample(1:nrow(x), 300),]
f3 <- x3[sample(1:nrow(x3), 250),];
y3 <- x3[!(x3$marker %in% f3$marker),];
y3$position <- as.numeric(y3$position);
y3$chromosome <- factor(y3$chromosome, levels=c(1:33, "X", "Y", "MT"));
y3 <- y3[order(y3$chromosome, y3$position), ];
write.table(x3, "cnfilt3a.in", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);
write.table(f3, "cnfilt3b.in", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);
write.table(y3, "cnfilt3.ans", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);
