#!/usr/bin/env python3
import sys
from pathlib import Path

if len(sys.argv) == 3:
	filename = Path(sys.argv[1])
	try:
		rate = float(sys.argv[2])
	except ValueError:
		print("USAGE: python3 gff3_genes_2_rmt.py [gff_file] [rate]")
		sys.exit(1)
else:
	print("USAGE: python3 gff3_genes_2_rmt.py [gff_file] [rate]")
	sys.exit(1)

rd = {}
try:
	with open(filename, "r") as hndl:
		for line in hndl.readlines():
			if line and not line.startswith("#"):
				values = line.strip().split("\t")
				if values[2].lower() == "gene":
					if values[0] not in rd.keys(
					):  # add lists under chromosome keys
						rd[values[0]] = []
					rd[values[0]].append([values[3],
							values[4]])  # add start stop values to lists

	with open(filename.with_suffix(".rmt"), "w") as hndl:
		hndl.write(f"std\nit None\nsn {rate}\n\n")
		chromosomes = list(rd.keys())
		for i in range(len(chromosomes)):
			hndl.write(f"chr {i+1} #{chromosomes[i]}\n")
			for interval in rd[chromosomes[i]]:
				hndl.write(f"{interval[0]}-{interval[1]} None\n")
except FileNotFoundError:
	print(f"ERROR: cannot find {filename}")
	sys.exit(1)
