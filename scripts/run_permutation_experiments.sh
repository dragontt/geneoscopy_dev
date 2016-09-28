#!/usr/bash
# Run LIMMA analysis on permuted expression data, sort top 60 DE genes, and analysis their results

# number of experiments
num_experiments=100

# permute gene expression data
original_expr_data=../data/20160128_project_1930/chipdata_geneset_x_valid_chips.txt
permuted_expr_data_dir=../data/20160128_project_1930/permutation_experiments
python permutate_expr_data.py -d ${original_expr_data} -n ${num_experiments} -o ${permuted_expr_data_dir}
echo "Generate perumted expression data."

# identify de genes
for i in $(seq 1 ${num_experiments}); do
	permuted_expr_data=../data/20160128_project_1930/permutation_experiments/permuted_expr_data_${i}.txt
	de_genes=../data/20160128_project_1930/permutation_experiments/60_limma_genes_${i}.txt
	Rscript de_analysis.r ${permuted_expr_data} ${de_genes}
	echo "Analyze DE genes on permutation experiments: ${i}..."
done

# analyze results
original_de_genes=../data/20160128_project_1930/60_limma_genes.txt
permuted_de_genes_dir=../data/20160128_project_1930/permutation_experiments/
analysis_results=../data/20160128_project_1930/permutation_experiments_results.txt
python permutate_de_gene_analysis.py -o ${original_de_genes} -d ${permuted_de_genes_dir} -n ${num_experiments} -a ${analysis_results}
echo "Analyze final results."
