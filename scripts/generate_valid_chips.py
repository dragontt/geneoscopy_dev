#!/usr/bin/python
import numpy as np

## new split 
proj_dir = '/Users/KANG/geneoscopy_dev/data/20170217_combined_projects_batch1-17/'
fn_expr = proj_dir + 'chipdata_header.txt'
fn_vc = proj_dir + 'valid_chips_x.txt'

f = open(fn_expr, 'r')
header = f.readline().strip().split("\t")
header = np.array(header[1:])
f.close()

total_c_te = 23
total_p_te = 20
total_n_te = 22

indx_c = []
indx_p = []
indx_n = []
for i in range(len(header)):
	if header[i].endswith('C'):
		indx_c.append(i)
	elif header[i].endswith('P'):
		indx_p.append(i)
	else:
		indx_n.append(i)

out = np.hstack((header[np.newaxis].T, 
	np.zeros(len(header), dtype=int)[np.newaxis].T, 
	np.zeros(len(header), dtype=int)[np.newaxis].T))
out[indx_c,1] = '2'
out[indx_p,1] = '1'

# out[np.random.choice(indx_c,size=total_c_te,replace=False),2] = '1'
# out[np.random.choice(indx_p,size=total_p_te,replace=False),2] = '1'
# out[np.random.choice(indx_n,size=total_n_te,replace=False),2] = '1'

indx_c = np.array(indx_c)
indx_p = np.array(indx_p)
indx_n = np.array(indx_n)
out[indx_c[np.linspace(0,len(indx_c)-1,num=total_c_te,dtype=int)],2] = '1'
out[indx_p[np.linspace(0,len(indx_p)-1,num=total_p_te,dtype=int)],2] = '1'
out[indx_n[np.linspace(0,len(indx_n)-1,num=total_n_te,dtype=int)],2] = '1'

np.savetxt(fn_vc, out, fmt="%s", delimiter="\t")



"""
## add split to existed split
proj_dir = '/Users/KANG/geneoscopy_dev/data/20170205_combined_projects_batch1-17/'
fn_in = proj_dir + '../20161115_combined_project_2283_abcdefghi/valid_chips.txt'
fn_expr = proj_dir + 'chipdata_geneset_x_valid_chips.txt'
fn_vc = proj_dir + 'valid_chips_x.txt'

f = open(fn_expr, 'r')
header = f.readline().strip().split("\t")
header = np.array(header[1:])
f.close()

vc_existed = np.loadtxt(fn_in, dtype=str, delimiter="\t")
out1 = []
c_te_existed = 0
p_te_existed = 0
n_te_existed = 0
for i in range(len(vc_existed)):
	[x,y,z] = vc_existed[i,]
	if x in header:
		out1.append(vc_existed[i,])
		if z == '1':
			if x.endswith('C'):
				c_te_existed += 1
			elif x.endswith('P'):
				p_te_existed += 1
			elif x.endswith('N'):
				n_te_existed += 1
out1 = np.array(out1)

total_c_tr = 76
total_p_tr = 62
total_n_tr = 66
total_c_te = 19
total_p_te = 16
total_n_te = 17

header_left = np.setdiff1d(header, out1[:,0])
indx_c = []
indx_p = []
indx_n = []
for i in range(len(header_left)):
	if header_left[i].endswith('C'):
		indx_c.append(i)
	elif header_left[i].endswith('P'):
		indx_p.append(i)
	else:
		indx_n.append(i)

out2 = np.hstack((header_left[np.newaxis].T, 
	np.zeros(len(header_left), dtype=int)[np.newaxis].T, 
	np.zeros(len(header_left), dtype=int)[np.newaxis].T))
out2[indx_c,1] = '2'
out2[indx_p,1] = '1'

# out2[np.random.choice(indx_c,size=total_c_te,replace=False),2] = '1'
# out2[np.random.choice(indx_p,size=total_p_te,replace=False),2] = '1'
# out2[np.random.choice(indx_n,size=total_n_te,replace=False),2] = '1'

indx_c = np.array(indx_c)
indx_p = np.array(indx_p)
indx_n = np.array(indx_n)
out2[indx_c[np.linspace(0,len(indx_c)-1,num=(total_c_te-c_te_existed),dtype=int)],2] = '1'
out2[indx_p[np.linspace(0,len(indx_p)-1,num=(total_p_te-p_te_existed),dtype=int)],2] = '1'
out2[indx_n[np.linspace(0,len(indx_n)-1,num=(total_n_te-n_te_existed),dtype=int)],2] = '1'

np.savetxt(fn_vc, np.vstack((out1, out2)), fmt="%s", delimiter="\t")
"""
