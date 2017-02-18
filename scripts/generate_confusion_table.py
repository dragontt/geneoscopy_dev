#!/usr/bin/python
import sys
import numpy as np

def calculate_confusion_matrix(label_te, label_pred, diseased_label):
	pred_pos = np.where(label_pred == diseased_label)[0]
	pred_neg = np.where(label_pred == "N")[0]
	te_pos = np.where(label_te == diseased_label)[0]
	te_neg = np.where(label_te == "N")[0]
	tps = len(np.intersect1d(pred_pos, te_pos)) 
	fps = len(np.intersect1d(pred_pos, te_neg)) 
	fns = len(np.intersect1d(pred_neg, te_pos))
	tns = len(np.intersect1d(pred_neg, te_neg))
	print tps, fps, fns, tns
	sens = tps/float(len(te_pos))
	spec = tns/float(len(te_neg))
	accu = (tps+tns)/float(len(label_te))
	return [sens, spec, accu]

file = sys.argv[1]
x = np.loadtxt(file, dtype=str, usecols=[1,2])
print calculate_confusion_matrix(x[:,0], x[:,1], "P")
print calculate_confusion_matrix(x[:,0], x[:,1], "C")

