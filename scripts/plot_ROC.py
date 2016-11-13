#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from sklearn import metrics

## file IO
dir_proj = "/Users/KANG/geneoscopy_dev/data/run_proj_abcdefgh/"
file_pred = dir_proj + "tmp/prediction_probs.txt"

f = open(file_pred, "r")
lines = f.readlines()

labels_all = []
labels_cancer = []
labels_polyp = []
labels_cvp = []
probs_all = []
probs_cancer = []
probs_polyp = []
probs_cvp = []

for i in range(len(lines)):
	line = lines[i].strip().split()
	pc, pn, pp = np.array(line[3:6], dtype=float)
	## positive vs negative recommendation 
	true_label = 0 if line[1] == "N" else 1
	labels_all.append(true_label)
	# probs_all.append(1-pn)
	probs_all.append(max(pc,pp) / (max(pc,pp)+pn))
	## cancer vs normal
	if line[1] == "N" or "C":
		true_label = 0 if line[1] == "N" else 1
		labels_cancer.append(true_label)
		probs_cancer.append(pc / (pc+pn))
	## polyp vs normal
	if line[1] == "N" or "P":
		true_label = 0 if line[1] == "N" else 1
		labels_polyp.append(true_label)
		probs_polyp.append(pp / (pp+pn))
	## cancer vs polyp
	if line[1] == "C" or "P":
		true_label = 0 if line[1] == "P" else 1
		labels_cvp.append(true_label)
		probs_cvp.append(pc / (pc+pp))

## plot multiple ROCs 
fpr_all, tpr_all, thresholds_all = metrics.roc_curve(labels_all, probs_all)
roc_auc_all = metrics.auc(fpr_all, tpr_all)
fpr_cancer, tpr_cancer, thresholds_cancer = metrics.roc_curve(labels_cancer, probs_cancer)
roc_auc_cancer = metrics.auc(fpr_cancer, tpr_cancer)
fpr_polyp, tpr_polyp, thresholds_polyp = metrics.roc_curve(labels_polyp, probs_polyp)
roc_auc_polyp = metrics.auc(fpr_polyp, tpr_polyp)
fpr_cvp, tpr_cvp, thresholds_cvp = metrics.roc_curve(labels_cvp, probs_cvp)
roc_auc_cvp = metrics.auc(fpr_cvp, tpr_cvp)

plt.figure()
plt.plot(fpr_all, tpr_all, color='r', linestyle='-', lw=4, 
	label='Pos/Neg recommendation (AUC = %0.2f)' % roc_auc_all)
plt.plot(fpr_cancer, tpr_cancer, color='g', linestyle='-', lw=2, 
	label='Cancer/Normal (AUC = %0.2f)' % roc_auc_cancer)
plt.plot(fpr_polyp, tpr_polyp, color='b', linestyle='-', lw=2, 
	label='Polyp/Normal (AUC = %0.2f)' % roc_auc_polyp)
plt.plot(fpr_cvp, tpr_cvp, color='m', linestyle='-', lw=2, 
	label='Cancer/Polyp (AUC = %0.2f)' % roc_auc_cvp)
plt.plot([0, 1], [0, 1], color='k', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(loc="lower right")
plt.show()

