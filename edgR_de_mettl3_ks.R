library("edgeR")

overlap_count_matrix <- "/Users/joshlevendis/rqc/28hpi_overlap_counts.txt"
x <- read.delim(overlap_count_matrix, header=TRUE, row.names="Geneid")
t1 <- x[,c(1, 2, 3, 4)]
groups <- factor(c(1, 1, 2, 2))
lib_sizes_28hpi <- c(2457554, 1362266, 1798829, 2299403)
y <- DGEList(counts=t1,group=groups, lib.sizes=lib_sizes_28hpi)

# fc_matrix <- "/Users/joshlevendis/Downloads/featureCountsStrandedOverlapMAPQFilter/8.3_featureCounts_headless"
# x <- read.delim(fc_matrix, header=TRUE, row.names="Geneid")
# t1 <- x[,c(6,7,8,9)]
# t2 <- x[,c(10,11,12,13)]
# t3 <- x[,c(14,15,16,17)]
# tall <- x[,c(6,7,8,9,10,11,12,13,14,15,16,17)]
# groups <- factor(c(1,1,2,2))
# groups_all <- factor(c(1,1,2,2,1,1,2,2,1,1,2,2))
# y <- DGEList(counts=t3,group=groups)
# 
keep <- filterByExpr(y)
#y <- y[keep,,keep.lib.sizes=FALSE]
#y <- normLibSizes(y)
design <- model.matrix(~groups)
y <- estimateDisp(y, design)


# exact test (edge v1)
qlf_exact <- exactTest(y)

# Negative binomial GLM (edge v2)
fit <- glmFit(y, design)
qlf_nbglm <- glmLRT(fit, coef=2)

# Quasi likelihood model (edge v3, v4)
fit <- glmQLFit(y, design)
qlf_qlm <- glmQLFTest(fit, coef=2)

qlf <- qlf_exact

# outfile <- "/Users/joshlevendis/Downloads/featureCountsStrandedOverlapMAPQFilter/edgeR_28hpi"
outfile <- "/Users/joshlevendis/rqc/edgeR_28hpi_overlaps"
tt = topTags(qlf,n=Inf)
df <- as.data.frame(tt)
colnames <- c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
plotMD(qlf)
write.table(df, file=outfile, sep="\t")

