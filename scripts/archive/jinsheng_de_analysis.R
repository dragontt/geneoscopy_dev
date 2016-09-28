library(limma)

data_table  <- read.delim(file="chipdata_table_xxx.txt", row.names = 1, header = T, sep = "\t", stringsAsFactors = F)

### your data table should have samples in coloumns and genes (probes) in rows, with annotations in extra columns

data_frame <- data_table[,c(xx:xy)] ## excluding annotation columns

Targets <- read.delim(file="target_xx_samples.txt", header = T, sep = "\t", stringsAsFactors = F)

f <- factor(Targets$Group, levels=c("A", "B", "C"))
design <- model.matrix(~0+f)
colnames(design) <- levels(f)

fit <- lmFit(data_frame, design)
contrast.matrix <- makeContrasts(AvsB=A-B, AvsC=A-C, BvsC=B-C, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

for (test in colnames(fit2$coefficients)[1:length(colnames(fit2$coefficients))]) { #{{
            results_name <- as.character(test)
            dir.create(results_name)
            o.dir <- setwd(results_name)
	    
            ### Do the differential analysis and multiple test correction ("adjust.method" defaut is BH)
            DIFF = topTable(fit2, sort = "none", n = Inf, confint = T, adjust.method="BH", coef = test)

            #Filter for significant results
            DIFF_SIG <-  DIFF[DIFF$P.Value <= 0.05, ]
            DIFF_FDR <-  DIFF[DIFF$adj.P.Val <= 0.05, ]

            #Write Data to a file
            OUT_DATA <- DIFF[order(DIFF$adj.P.Val), ]
            OUT_DATA = format(OUT_DATA, round = 5)
            write.table(cbind(OUT_DATA, data_table[match(rownames(OUT_DATA), rownames(data_table)), , drop=FALSE]), file=paste(results_name,"Differential_Expression.txt", sep = "_"), col.names=NA, row.names=TRUE, sep = "\t", quote = FALSE)

	    setwd(o.dir)

} #}}
