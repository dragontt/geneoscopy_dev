#!/usr/bin/python
import sys
import argparse
import numpy as np

def parse_args(argv):
	parser = argparse.ArgumentParser(description='Perform quality control.')
	parser.add_argument('-i', '--input_expr', dest='input_expr')
	parser.add_argument('-l', '--transcript_cluster_list', dest='transcript_cluster_list')
	parser.add_argument('-c', '--column_number', dest='column_number', type=int)
	parser.add_argument('-o', '--output_expr', dest='output_expr')
	parsed = parser.parse_args(argv[1:])
	return parsed

def main(argv):
	parsed = parse_args(argv)

	f = open(parsed.transcript_cluster_list, "r")
	lines = f.readlines()
	tcs = set()
	for line in lines:
		tmp_arr = line.strip().split("\t")
		if len(tmp_arr) > parsed.column_number:
			for tmp_tc in tmp_arr[parsed.column_number].split(","):
				tcs.add(tmp_tc)
	f.close()
	tcs = list(tcs)
	print "Transript cluster count:", len(tcs)

	expr = np.loadtxt(parsed.input_expr, dtype=str, delimiter="\t")
	tc_indx = [0]
	for i in range(1,len(expr)):
		if expr[i,0] in tcs:
			tc_indx.append(i)

	np.savetxt(parsed.output_expr, expr[tc_indx,:], fmt="%s", delimiter="\t")

if __name__ == "__main__":
    main(sys.argv)
