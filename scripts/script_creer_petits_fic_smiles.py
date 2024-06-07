#!/usr/bin/env python
import csv
import sys
import os

file_data = sys.argv[1]
with open(file_data, 'r') as f :
	csvreader = csv.reader(f)
	i = 0
	for row in csvreader :
		if row[1] != "name":
			petit_smile = row[2]
			name = "_".join(row[1].split(" "))
			with open("data/CHIMISTE/smi_files/%s.smi"%name, 'w') as f :
				f.write("%s %s\n"%(petit_smile, name))
			i += 1
