## Example:
## Rscript remove_batch_effect.r chipdata_rma.expression_console.txt sample_batch.txt chipdata_rma.expression_console.batcheff_removed.txt

require(sva)

## parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
file_expr <- args[1]
file_batch <- args[2]
file_out <- args[3]

## load and prepare data
expr <- read.table(file_expr, sep="\t", header=T, row.names=1, check.names=F)
batch <- read.table(file_batch, sep="\t", header=T, check.names=F)$Batch
mod <- model.matrix(~1, data=as.data.frame(batch))

## remove batch using ComBat
expr_out <- ComBat(dat=expr, batch=batch, mod=mod)

## output expression data
write.table(expr_out, file_out, sep="\t", quote=F, col.names=NA)