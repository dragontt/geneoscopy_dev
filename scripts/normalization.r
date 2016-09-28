### Usage ###
# Rscript rma normalization.r <dir_.cel_files> <normalized_expr_file>

## download and load libraries
source("http://bioconductor.org/biocLite.R")
if (!require("oligo")) try(biocLite("oligo"));
if (!require("limma")) try(biocLite("limma"));
if (!require("gcrma")) try(biocLite("gcrma"));
if (!require("hta20sttranscriptcluster.db")) try(biocLite("hta20sttranscriptcluster.db"));

## File IO
args <- commandArgs(trailingOnly = TRUE)
norm_method <- args[1]
dir_raw <- args[2]
file_expr <- args[3]
# file_annot <- args[4]

if (norm_method == "rma") {
	library(limma)
	library(oligo)
	library(hta20sttranscriptcluster.db)

	## rma normalization  
	celFiles <- list.celfiles(dir_raw, full.names=T)
	affyRaw <- read.celfiles(celFiles)
	eset <- rma(affyRaw, target="core")
	# eset <- rma(affyRaw, target="probeset")
	write.exprs(eset,file=file_rma_expr)

	## add gene annotation
	# my_frame <- data.frame(exprs(eset))
	# Annot <- data.frame(ACCNUM=sapply(contents(hta20sttranscriptclusterACCNUM),paste,collapse=","), SYMBOL=sapply(contents(hta20sttranscriptclusterSYMBOL),paste,collaps=","), DESC=sapply(contents(hta20sttranscriptclusterGENENAME),paste,collapse=","))
	# all <- merge(Annot, my_frame, by.x=0, by.y=0, all=T)
	# write.table(all,file=file_annot,sep="\t",quote=FALSE)
} else if (norm_method == "gcrma") {
	library(oligo)
	library(gcrma)
	## gcrma normalization
	celFiles <- list.celfiles(dir_raw, full.names=T)
	# affyRaw <- read.celfiles(celFiles)
	affyRaw <- ReadAffy()
	eset <- gcrma(affyRaw)
	write.exprs(eset,file=file_gcrma_expr)
} else {
	cat("No normalization method!", "\n")
}
