library(rgl)
## Example:
## Rscript pca_analysis.r /Users/KANG/geneoscopy_dev/data/run_proj_a_b_c_d/training/chipdata.txt /Users/KANG/geneoscopy_dev/data/run_proj_a_b_c_d/training/top_de_genes.txt /Users/KANG/geneoscopy_dev/data/run_proj_a_b_c_d/training/valid_chips.txt

## Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
file_expr_data <- args[1]
file_gene_list <- args[2]
file_sample_labels <- args[3]

## Run principal component analysis
expr_data <- read.table(file_expr_data, sep='\t', row.names=1, header=TRUE, fill=TRUE)
expr_data[is.na(expr_data)] <- 0

## Filter based on top genes
top_genes <- read.table(pipe(paste('cut -f1', file_gene_list)), skip=1)
row_idx <- c()
for (i in 1 : dim(top_genes)[1]) {
	matched_idx <- match(top_genes[i,], row.names(expr_data))
	row_idx <- append(row_idx, matched_idx)
}
expr_data <- expr_data[row_idx,]

## pca analysis
expr_pca <- prcomp(t(expr_data))

## Plot 2d pca
# biplot(expr_pca,var.axes=FALSE)

## Plot 3d pca
plot3d(expr_pca$x,xlab="PC1",ylab="PC2",zlab="PC3",type="h")
# spheres3d(expr_pca$expr_data, radius=3, col=rainbow(length(expr_pca$ex[,1])))
# colors <- list(n='chartreuse3', p='goldenrod1', c='firebrick1')
# label_colors <- c(rep(colors$n,5), colors$p, rep(colors$c,4), colors$p, colors$p)
colors <- list('0'='chartreuse3', '1'='firebrick1')
labels <- read.table(file_sample_labels, sep="\t", colClasses=c("NULL","character","NULL"))
label_colors <- c()
for (i in 1:dim(labels)[1]){
	label <- labels[i,1]
	label_colors <- c(label_colors, colors[[label]])
}
spheres3d(expr_pca$x, radius=.25, col=label_colors)
grid3d(side="z", at=list(z=0))
# text3d(expr_pca$x, text=rownames(expr_pca$x), adj=1.3)

## Pause after PCA plot
