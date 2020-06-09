#!/usr/bin/env python3
import sys

if len(sys.argv) > 1:
	filename = sys.argv[1]
else:
	print("ERROR: Please specify the GFF3 file.")

rd = {}
with open(filename, "r") as hndl:
	for line in hndl.readlines():
		if line and not line.startswith("#"):
			values = line.strip().split("\t")
			if values[2].lower() == "gene":
				if values[0] not in rd.keys(): # add lists under chromosome keys
					rd[values[0]] = []
				rd[values[0]].append([values[3], values[4]]) # add start stop values to lists

with open(".".join(filename.split(".")[:-1])+".rmt", "w") as hndl:
	hndl.write("std\nit None\nsn 0.01\n\n")
	chromosomes = list(rd.keys())
	for i in range(len(chromosomes)):
		hndl.write(f"chr {i+1} #{chromosomes[i]}\n")
		for interval in rd[chromosomes[i]]:
			hndl.write(f"{interval[0]}-{interval[1]} None\n")
