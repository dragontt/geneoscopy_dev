library(rgl)

# run principal component analysis
expr_data <- read.table('/Users/KANG/geneoscopy_dev/data/20160808_project_2283/normal_vs_polyp_vs_crc/chipdata_geneset_x_valid_chips.txt', sep='\t', row.names=1, header=TRUE, fill=TRUE)
expr_data[is.na(expr_data)] <- 0

# filter based on top genes
num_genes <- 59
top_genes <- read.table(pipe('cut -f1 /Users/KANG/geneoscopy_dev/data/20160808_project_2283/normal_vs_polyp_vs_crc/top_59_limma_genes_pval_0.001.txt'), skip=1)
# num_genes <- 1000
# top_genes <- read.table(pipe('cut -f1 /Users/KANG/geneoscopy_dev/data/20160808_project_2283/normal_vs_polyp_vs_crc/top_67529_limma_genes.txt'), skip=1)

row_idx <- c()
for (i in 1 : max(dim(top_genes)[1], num_genes)) {
	matched_idx <- match(top_genes[i,], row.names(expr_data))
	row_idx <- append(row_idx, matched_idx)
}
expr_data <- expr_data[row_idx,]

# pca analysis
expr_pca <- prcomp(t(expr_data))

# # plot 2d pca
# biplot(expr_pca,var.axes=FALSE)

# plot 3d pca
plot3d(expr_pca$x,xlab="PC1",ylab="PC2",zlab="PC3",type="h")
# spheres3d(expr_pca$expr_data, radius=3,col=rainbow(length(expr_pca$ex[,1])))
colors <- list(n='chartreuse3', p='goldenrod1', c='firebrick1')
labels <- c(rep(colors$n,5), colors$p, rep(colors$c,4), colors$p, colors$p)
spheres3d(expr_pca$x, radius=.25,col=labels)
grid3d(side="z", at=list(z=0))
text3d(expr_pca$x, text=rownames(expr_pca$x), adj=1.3)

