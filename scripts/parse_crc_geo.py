#/usr/bin/python
import numpy

path_proj = "/Users/KANG/geneoscopy/"
path_out = path_proj + "output/crc_geo_data/GDS4382/"
file_in = path_proj + "resources/crc_geo_data/GDS4382_full.soft"

lines = open(file_in, "r").readlines()

feature = []
mtr = []

pid_crc = lines[26].split('=')[1].strip().split(',')
pid_norm = lines[31].split('=')[1].strip().split(',')
pid = lines[181].split()[2:36]
label = [None] * 34
for i in range(len(pid)):
	if pid[i] in pid_crc:
		label[i] = 1
	elif pid[i] in pid_norm:
		label[i] = 0
mtr.append(label)

# for i in range(182, len(lines)):
for i in range(182,len(lines)-1):
	line = lines[i].split()
	feature.append(line[0:2])
	mtr.append(line[2:36])

pid = numpy.transpose(numpy.array(pid))
feature = numpy.array(feature)
mtr = numpy.transpose(numpy.array(mtr))

numpy.savetxt(path_out + 'patient_id.txt', pid, fmt='%s', delimiter='\t', newline='\n')
numpy.savetxt(path_out + 'feature_id.txt', feature, fmt='%s', delimiter='\t', newline='\n')
numpy.savetxt(path_out + 'expression.mtr', mtr, fmt='%s', delimiter='\t', newline='\n')