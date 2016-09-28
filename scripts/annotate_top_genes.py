import sys
import argparse
import numpy as np
import json

def parseArgs(argv):
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--input_table', dest='input_table')
    parser.add_argument('-j', '--json_dict', dest='json_dict')
    parser.add_argument('-o', '--output_table', dest='output_table')
    parsed = parser.parse_args(argv[1:])
    return parsed

def concatenateList(data):
	if len(data) > 0:
		return ",".join(data)
	else:
		return "---"

def chooseForDisplay(gene_sym, gene_acc, rna_acc, probeset):
	if len(gene_sym) > 0:
		return gene_sym[0]
	elif len(gene_acc) > 0:
		return gene_acc[0]
	elif len(rna_acc) > 0:
		return rna_acc[0]
	else:
		return probeset

def main(argv):
	## Load data
	parsed = parseArgs(argv)
	with open(parsed.json_dict) as json_data:
		probeset_dict = json.load(json_data)
	top_genes = np.loadtxt(parsed.input_table, dtype=str, skiprows=1, usecols=[0,1,4])

	## Annotate top genes table
	writer = open(parsed.output_table, "w")
	writer.write("probeset\tgene_symbol\tgene_accession\trna_accession\tfor_display\tlogFC\tp_value\n")
	for i in range(len(top_genes)):
		probeset = top_genes[i,0]
		probeset_attributes = probeset_dict[probeset]
		gene_sym = probeset_attributes["gene_symbol"]
		gene_acc = probeset_attributes["gene_acc"]
		rna_acc = probeset_attributes["rna_acc"]
		for_display = chooseForDisplay(gene_sym, gene_acc, rna_acc, probeset)
		gene_sym = concatenateList(gene_sym)
		gene_acc = concatenateList(gene_acc)
		rna_acc = concatenateList(rna_acc)
		writer.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" 
			% (probeset, gene_sym, gene_acc, rna_acc, for_display, top_genes[i,1], top_genes[i,2]))
	writer.close()

if __name__ == "__main__":
    main(sys.argv)
