#!/usr/bin/python
import sys
import numpy as np
from random import shuffle

file_data = sys.argv[1]
file_rnd_genes = sys.argv[2]

gids = np.loadtxt(file_data, dtype=str, skiprows=1, usecols=[0])
indx = range(len(gids))
shuffle(indx)
np.savetxt(file_rnd_genes, gids[indx[:200]], fmt='%s')
