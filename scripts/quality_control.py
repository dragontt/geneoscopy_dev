#!/usr/bin/python
# example:
# python quality_control.py -n 15 -d ../data/20160128_project_1930/chipdata_HTA2.0_15chips_Barnell_Project1930.RMA-GENE-FULL.TXT -q ../data/20160128_project_1930/QC_Table_HTA2.0_15chips_Barnell_Project1930.txt -s ../data/20160128_project_1930/sample_sheet_project_1930.txt -v ../data/20160128_project_1930/valid_chips.txt -o ../data/20160128_project_1930/chipdata_geneset_x_valid_chips.txt

import sys
import argparse
import numpy
from random import random

label_dict = {"N_vs_C":{"N":0, "C":1}, "N_vs_P_vs_C":{"N":0, "P":1, "C":2}}

def parse_args(argv):
	parser = argparse.ArgumentParser(description='Perform quality control.')
	parser.add_argument('-n', '--num_samples', dest='num_samples', type=int)
	parser.add_argument('-d', '--data', dest='data')
	parser.add_argument('-q', '--qc_table', dest='qc_table')
	parser.add_argument('-s', '--sample_sheet', dest='sample_sheet')
	parser.add_argument('-g', '--group', dest='group')
	parser.add_argument('-v', '--valid_chips', dest='valid_chips')
	parser.add_argument('-o', '--valid_data', dest='valid_data')
	parsed = parser.parse_args(argv[1:])
	return parsed


def filter_rle_mean(qc_table, rle_mean_threshold):
	valid_chips = []
	for i in range(len(qc_table)):
		sample = qc_table[i,0].split('.')[0]
		rle_mean = float(qc_table[i,4])
		if rle_mean <= rle_mean_threshold:
			valid_chips.append(sample)
		valid_chips.append(sample)
	return valid_chips


def filter_pos_vs_neg_auc(qc_table, auc_threshold):
	valid_chips = []
	for i in range(len(qc_table)):
		sample = qc_table[i,0].split('.')[0]
		auc = float(qc_table[i,5])
		if auc >= auc_threshold:
			valid_chips.append(sample)
	return valid_chips


def compute_group_difference():
	pass


def annotate_samples(valid_chips, sample_sheet):
	valid_chips_out = []
	for i in range(len(sample_sheet)):
		chip_id = sample_sheet[i,2]
		sample_id = sample_sheet[i,3]
		if len(sample_id.split(".")) < 2:
			sys.exit("WARNING: ChIP ID "+ chip_id +" does not have proper label.") 
		if chip_id in valid_chips:
			valid_chip_idx = valid_chips.index(chip_id)
			valid_chip = valid_chips[valid_chip_idx]
			valid_chips_out.append(".".join([valid_chip, sample_id.split(".")[1]]))
	return valid_chips_out


def filter_data(valid_chips, data):
	valid_col_idx = [0]
	valid_row_idx = []

	data_header = [data[0,0]]
	for i in range(1, len(data[0])):
		chip_id = data[0,i].split('.')[0]
		for valid_chip in valid_chips:
			if chip_id == valid_chip.split('.')[0]:
				valid_col_idx.append(i)
				data_header.append(valid_chip)
	data_header = numpy.array(data_header, dtype=str)

	for i in range(1, len(data[:,0])):
		if data[i,0].startswith('TC'):
			valid_row_idx.append(i)

	new_data = data[valid_row_idx]
	new_data = new_data[:,valid_col_idx]
	new_data = numpy.vstack((data_header, new_data))
	return new_data


def write_output(fn_valid_chips, fn_data, valid_chips, data, group):
	valid_chips = numpy.array(valid_chips, dtype=str)[numpy.newaxis].T
	labels = -1*numpy.ones(valid_chips.shape, dtype=int)
	test_flag = -1*numpy.ones(valid_chips.shape, dtype=int)
	for i in range(valid_chips.shape[0]):
		tmp = valid_chips[i,0].split(".")[1]
		labels[i] = label_dict[group][tmp]
		test_flag[i] = 1 if random() > 0.8 else 0
	# print "# of samples for testing:", numpy.count_nonzero(test_flag)
	# numpy.savetxt(fn_valid_chips, numpy.hstack((valid_chips, labels, test_flag)), delimiter='\t', fmt='%s')
	numpy.savetxt(fn_data, data, delimiter='\t', fmt='%s')


def main(argv):
	parsed = parse_args(argv)
	num_samples = parsed.num_samples

	# load files
	data = numpy.loadtxt(parsed.data, dtype=str, delimiter='\t', usecols=range(num_samples+1))
	qc_table = numpy.loadtxt(parsed.qc_table, dtype=str, delimiter='\t', skiprows=1)
	sample_sheet = numpy.loadtxt(parsed.sample_sheet, dtype=str, delimiter='\t', skiprows=1)
	
	# filter out samples with rle mean > 0.23~0.25
	# valid_chips = filter_rle_mean(qc_table, .25)

	# group samples by pos vs neg auc
	valid_chips = filter_pos_vs_neg_auc(qc_table, .75)

	# check if auc mean of diseased vs mean of normal is sig different
	# pass

	# append normal, polyps or CRC labels to the normal sample
	# and randomize 80% for training and 205 for testing
	valid_chips = annotate_samples(valid_chips, sample_sheet)

	# filter data using valid samples, and only use chipset with TC prefix
	filtered_data = filter_data(valid_chips, data)

	# write valid samples and data
	write_output(parsed.valid_chips, parsed.valid_data, valid_chips, filtered_data, parsed.group)

if __name__ == "__main__":
    main(sys.argv)
