#!/usr/bin/python
import numpy as np

threshold = 4

dir_proj = "/Users/KANG/geneoscopy_dev/data/run_proj_batch1-17_4/"
filename_stools = dir_proj + "samples_w_gram.txt"
filename_vc = dir_proj + "valid_chips.txt"
filename_vc_out = dir_proj + "valid_chips.high_stool.txt"

# labels = {}
high_stool = []
f = open(filename_stools, "r")
lines = f.readlines()
for i in range(len(lines)):
	line = lines[i].strip().split("\t")
	if line[0].startswith("2"):
		if float(line[1]) >= threshold:
			high_stool.append(line[0])
		# if line[2] == "Normal":
		# 	labels[line[0]] = '0'
		# elif line[2] == "Polyp":
		# 	labels[line[0]] = '1'
		# elif line[2] == "Cancer":
		# 	labels[line[0]] = '2'
f.close()

vc = np.loadtxt(filename_vc, dtype=str)
indx = []
for i in range(len(vc)):
	chipid = vc[i,0].split(".")[0]
	if chipid in high_stool:
		indx.append(i)
	# if labels[chipid] != vc[i,1]:
	# 	print labels[chipid], vc[i,0], vc[i,1]

print "N", len(np.where(vc[indx,1] == '0')[0])
print "P", len(np.where(vc[indx,1] == '1')[0])
print "C", len(np.where(vc[indx,1] == '2')[0])

np.savetxt(filename_vc_out, vc[indx,], fmt="%s", delimiter="\t")
