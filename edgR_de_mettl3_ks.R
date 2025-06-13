library("edgeR")
fc_matrix <- "/Users/joshlevendis/Downloads/featureCountsStranded/8.3_featureCounts_headless"

x <- read.delim(fc_matrix, header=TRUE, row.names="Geneid")

t1 <- x[,c(6,7,8,9)]
t2 <- x[,c(10,11,12,13)]
t3 <- x[,c(14,15,16,17)]
groups <- factor(c(1,1,2,2))
y <- DGEList(counts=t3,group=groups)


keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- normLibSizes(y)
design <- model.matrix(~groups)
fit <- glmQLFit(y, design, robust=TRUE)

plotQLDisp(fit)

qlf <- glmQLFTest(fit, coef=2)

outfile <- "/Users/joshlevendis/Downloads/featureCountsStranded/edgeR_36hpi"
tt = topTags(qlf,n=Inf)
df <- as.data.frame(tt)
colnames <- c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
#write.table(df, file=outfile, sep="\t")

