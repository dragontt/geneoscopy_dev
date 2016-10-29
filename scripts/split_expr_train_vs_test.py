#!/usr/bin/python
import sys
import argparse
import numpy as np

def parse_args(argv):
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-v', '--valid_chips', dest='valid_chips')
	parser.add_argument('-i', '--input_expr', dest='input_expr')
	parser.add_argument('-tr0', '--train_expr', dest='train_expr')
	parser.add_argument('-tr1', '--train_valid_chips', dest='train_valid_chips')
	parser.add_argument('-te0', '--test_expr', dest='test_expr')
	parser.add_argument('-te1', '--test_valid_chips', dest='test_valid_chips')
	parsed = parser.parse_args(argv[1:])
	return parsed

def main(argv):
	parsed = parse_args(argv)
	
	valid_chips = np.loadtxt(parsed.valid_chips, dtype=str)
	indx_tr = np.where(valid_chips[:,2] == "0")[0]
	indx_te = np.where(valid_chips[:,2] == "1")[0]

	expr = np.loadtxt(parsed.input_expr, dtype=str, delimiter="\t")
	rownames = expr[:,0][np.newaxis].T
	expr = expr[:,1:]

	np.savetxt(parsed.train_expr, np.hstack((rownames, expr[:,indx_tr])), fmt="%s", delimiter="\t")
	if parsed.train_valid_chips != None:
		np.savetxt(parsed.train_valid_chips, valid_chips[indx_tr,:], fmt="%s", delimiter="\t")
	np.savetxt(parsed.test_expr, np.hstack((rownames, expr[:,indx_te])), fmt="%s", delimiter="\t")
	if parsed.test_valid_chips != None:
		np.savetxt(parsed.test_valid_chips, valid_chips[indx_te,:], fmt="%s", delimiter="\t")

if __name__ == "__main__":
    main(sys.argv)