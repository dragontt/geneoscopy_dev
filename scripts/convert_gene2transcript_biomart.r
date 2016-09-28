#!/usr/bin/R
library('biomaRt')

args <- commandArgs(trailingOnly = TRUE)
filename_input <- args[1]
use_col <- as.integer(args[2])
filename_output <- args[3]

input <- as.matrix(read.table(filename_input, sep="\t"))
genes <- unique(input[,use_col])

output <- c()
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "hgnc_symbol", 
	attributes=c("ensembl_gene_id","ensembl_transcript_id","hgnc_symbol"), 
	values=genes, mart=mart)

for (i in 1:length(genes)) {
	gene <- genes[i]
	ind = which(G_list[,3] == gene)
	output <- rbind(output, c(gene, 
		paste(unique(G_list[ind,1]), collapse=","), 
		paste(G_list[ind,2], collapse=","), 
		paste(input[which(input[,use_col] == gene),1], collapse=",")))
}

write.table(output, file=filename_output, row.names=F, col.names=F, quote=F, sep="\t")
