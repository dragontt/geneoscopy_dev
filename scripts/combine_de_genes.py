#!/usr/bin/python
import sys
import numpy as np

file_constr = sys.argv[1]
file_all = sys.argv[2]
pval_constr = float(sys.argv[3])
pval_all = float(sys.argv[4])
file_out = sys.argv[5]

constr_de = np.loadtxt(file_constr, dtype=str, skiprows=1, usecols=[0,4])
all_de = np.loadtxt(file_all, dtype=str, skiprows=1, usecols=[0,4])

tcs_c = constr_de[np.array(constr_de[:,1], dtype=float) <= pval_constr, 0]
tcs_a = all_de[np.array(all_de[:,1], dtype=float) <= pval_all, 0]
out = np.append("TC", np.union1d(tcs_c, tcs_a))

np.savetxt(file_out, out, fmt="%s")
