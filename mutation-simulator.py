#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Mutation-Simulator Version min-mut-len-2
# Copyright (C) 2019 Marius Kühl

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from tqdm import tqdm
from numpy.random import choice, random, randint
from numpy import arange, setdiff1d, array
from blist import blist
from pyfaidx import Fasta, Faidx
import random as rnd
import sys
import hashlib
import datetime
import argparse
import timeit
import itertools


def main():
	"""The main function."""
	start = timeit.default_timer()  # runtime
	filename, outfile_basename, rmt, mut_rates, mut_lengs, mut_block, assembly_name, species_name, sample_name, it_rate, ignore_warnings, titv = utilise_sysargs()  # checks the commandline arguments
	if not filename:
		return False
	fasta = read_fasta(filename)  # returns pyfaidx's fasta
	if not fasta:
		return False
	fai = read_fai(filename)
	if not fai:
		return False
	if rmt:
		rmt, mut_flag, it_flag, species_name, assembly_name, sample_name, mut_block, titv = load_rmt(rmt, fai, filename, ignore_warnings)
		if not rmt:
			return False
	else:
		rmt, mut_flag, it_flag = None, None, None
	if mut_rates or mut_flag:
		if not mutator(fasta, fai, mut_rates, mut_lengs, mut_block, rmt, outfile_basename, filename, assembly_name, species_name, sample_name, titv):
			return False
	if it_rate or it_flag:
		records_it, it_changes = interchromosomal_transloc(fasta, it_rate, rmt)
		if records_it and it_changes:
			if not subdicts_empty(it_changes):
				mode ="w"
				chromosomes = list(fasta.keys())
				for i in range(len(chromosomes)):
					write_fasta(outfile_basename+".fa", records_it[i], fasta[chromosomes[i]].long_name, fai.index[chromosomes[i]].lenc, mode)
					mode = "a"
				save_it_bedpe(outfile_basename+".bedpe", fasta, it_changes)
			else:
				print("ERROR: No interchromosomal translocations could be generated. Rate is probably too low.")
	stop = timeit.default_timer()
	runtime = round(stop - start, 4)
	sec = datetime.timedelta(seconds=runtime)
	d = datetime.datetime(1, 1, 1) + sec
	print(f"\nMutation-Simulator finished in: {runtime}s -> {d.hour}h {d.minute}m {d.second}s")
	return True


def read_fai(filename):
	"""Creates and reads a Fasta index file and returns content as a pyfaidx.Faidx"""
	try:
		return Faidx(filename, one_based_attributes=False)
	except(FileNotFoundError, IOError):
		print(f"ERROR: {filename} cannot be found")
		return False
	except ValueError:
		print(f"ERROR: Fasta file {filename} includes duplicate header")
		return False


def read_fasta(filename):
	"""Reads a Fasta file and returns its content as a pyfaidx.Fasta."""
	try:
		return Fasta(filename, one_based_attributes=False, as_raw=True, sequence_always_upper=True, read_ahead=10000)
	except(FileNotFoundError, IOError):
		print(f"ERROR: {filename} cannot be found")
		return False
	except ValueError:
		print(f"ERROR: Fasta file {filename} includes duplicate header")
		return False


def mutator(fasta, fai, mut_rates, mut_lengs, mut_block, rmt, outfile_basename, ref_filename, assembly_name, species_name, sample_name, titv):
	"""
	Generates and implements mutations on a given Chromosome with rate or RMT information.
	If rmt is given, mut_rates, mut_lengs and mut_block will be ignored. RMT must be set to False or None if it
	should not be used.

	Parameters:
		fasta (pyfaidx.Fasta): Chromosome information from any Fasta index file.
		fai (pyfaidx.Faidx): Fasta index file content.
		mut_rates (dict): All rates for all mutation types with 2 letter acronymes.
		mut_lengs (dict): All mutation lengths settings for any type.
		mut_block (dict): The range of the blocked area after any mutation for all types of mutations.
		rmt (dict): The rmt information from load_rmt(). Can be None.
		outfile_basename (str): Basename for the Fasta and VCF output file.
		assembly_name (str): Name of the Assembly.
		species_name (str): Name of the species.
		sample_name (str): Name of the sample.
		titv (float): Transition / Transversion ratio.
	"""
	chromosomes = list(fasta.keys())
	has_written_flag = False
	pbar = tqdm(total=len(chromosomes), desc="Mutating Sequences", position=0)
	for i in range(len(chromosomes)):
		no_tl_regions = []
		if rmt:
			mut_list = blist()
			translocations = blist()
			blocked_positions = set()
			for n in range(0, len(rmt[i][0])):
				if rmt[i][0][n][1]:
					range_mut_list = get_mutations(rmt[i][0][n][0][0], rmt[i][0][n][0][1], rmt[i][0][n][1], rmt[i][0][n][2], mut_block)
					if range_mut_list[0]:
						if range_mut_list[0] == -1:
							no_tl_regions.append(range(rmt[i][0][n][0][0], rmt[i][0][n][0][1] + 1))
						else:
							mut_list = mut_list + range_mut_list[0]
							translocations = translocations + range_mut_list[1]
							blocked_positions.update(range_mut_list[2])
					else:
						return False
					if not rmt[i][0][n][3]:
						no_tl_regions.append(range(rmt[i][0][n][0][0], rmt[i][0][n][0][1] + 1))
				else:
					no_tl_regions.append(range(rmt[i][0][n][0][0], rmt[i][0][n][0][1] + 1))
		else:
			mut_list, translocations, blocked_positions = get_mutations(0, len(fasta[chromosomes[i]])-1, mut_rates, mut_lengs, mut_block)
			if not mut_list:
				return False
		if mut_list != -1:
			if not sublists_empty(mut_list):
				if translocations:
					mut_list = mut_list + get_trans_inserts(translocations, len(fasta[chromosomes[i]]), blocked_positions, no_tl_regions)
				del translocations[:]
				del translocations
				del blocked_positions
				mut_list.sort(key=lambda x: x[1], reverse=True)
				record, mut_list = mutate(fasta[chromosomes[i]], mut_list, titv)
				if not has_written_flag:
					mode = "w"
				else: mode = "a"
				write_fasta(outfile_basename + ".fa", record, fasta[chromosomes[i]].long_name, fai.index[chromosomes[i]].lenc, mode)
				mut_list.sort(key=lambda x: x[1])
				save_mutations_vcf(outfile_basename + ".vcf", ref_filename, fasta, chromosomes[i], mut_list, assembly_name, species_name, sample_name, mode)
				has_written_flag = True
				del record
				del mut_list[:]
				del mut_list
			else:
				chr_has_rates = False
				for y in range(len(rmt[i][0])):
					if rmt[i][0][y][1]:
						chr_has_rates = True
				if chr_has_rates:
					pbar.write(f"WARNING: No mutations could be generated on sequence {i + 1}. (Mutation rate too low)")
		else:
			pbar.write(f"WARNING: No mutations could be generated on sequence {i + 1}. (Mutation rate too low)")
		pbar.update(1)
	pbar.close()
	return True


def utilise_sysargs():
	"""
	Uses argparse to utilise system arguments and provides default values.

	Returns:
		filename (str): Name of the input Fasta file.
		outfile_basename (str): Basename for the Fasta and VCF output file.
		rmt_file (str): Name of the RMT file. Can be None.
		mut_rates (dict): All rates for all mutation types with 2 letter acronymes.
		mut_lengs (dict): All mutation lengths settings for any mutation type.
		mut_block (dict): The range of the blocked area after any mutation for all types of mutations.
		assembly_name (str): Name of the assembly.
		species_name (str): Name of the species.
		sample_name (str): Name of the sample.
		it_rate(float): Rate at which interchromosomal translocations occure.
		ignore_warnings (bool): Ignores RMT warnings.
		titv (float): Transition / Transversion ratio.
	"""
	parser = argparse.ArgumentParser("See https://github.com/mkpython3/Mutation-Simulator/blob/master/README.md for more information about this program.")
	parser.add_argument("file", help="Path of the reference Fasta file")
	parser.add_argument("-o", "--output", help="Path and basename for the output files (without file extension)", default="")
	subparsers = parser.add_subparsers(help="Generate mutations or interchromosomal translocations via rmt or arguments")
	parser_rmt = subparsers.add_parser("rmt", help="Use random mutation table instead of arguments")
	parser_rmt.add_argument("rmtfile", help="Path to the rmt file")
	parser_rmt.add_argument("--ignore-warnings", help="Ignores RMT warnings.", action='store_true', default=False)
	parser_args = subparsers.add_parser("args", help="Use commandline arguments for mutations instead of rmt")
	parser_args.add_argument("-sn", "--snp", help="SNP rate. Default = 0", type=float, default=0)
	parser_args.add_argument("-snb", "--snpblock", help="Amount of bases blocked after SNP. Default = 1", type=int, default=1)
	parser_args.add_argument("-titv", "--transitionstransversions", help="Ratio of transitions:transversions likelihood. Default = 1", type=float, default=1)
	parser_args.add_argument("-in", "--insert", help="Insert rate. Default = 0", type=float, default=0)
	parser_args.add_argument("-inmin", "--insertminlength", help="Minimum length of inserts. Default = 1", type=int, default=1)
	parser_args.add_argument("-inmax", "--insertmaxlength", help="Maximum length of inserts. Default = 2", type=int, default=2)
	parser_args.add_argument("-inb", "--insertblock", help="Amount of bases blocked after insert. Default = 1", type=int, default=1)
	parser_args.add_argument("-de", "--deletion", help="Deletion rate. Default = 0", type=float, default=0)
	parser_args.add_argument("-demin", "--deletionminlength", help="Minimum length of deletions. Default = 1", type=int, default=1)
	parser_args.add_argument("-demax", "--deletionmaxlength", help="Maximum length of deletions. Default = 2", type=int, default=2)
	parser_args.add_argument("-deb", "--deletionblock", help="Amount of bases blocked after deletion. Default = 1", type=int, default=1)
	parser_args.add_argument("-iv", "--inversion", help="Inversion rate. Default = 0", type=float, default=0)
	parser_args.add_argument("-ivmin", "--inversionminlength", help="Minimum length of inversion. Default = 2", type=int, default=2)
	parser_args.add_argument("-ivmax", "--inversionmaxlength", help="Maximum length of inversion. Default = 3", type=int, default=3)
	parser_args.add_argument("-ivb", "--inversionblock", help="Amount of bases blocked after inversion. Default = 1", type=int, default=1)
	parser_args.add_argument("-du", "--duplication", help="Duplication rate. Default = 0", type=float, default=0)
	parser_args.add_argument("-dumin", "--duplicationminlength", help="Minimum length of duplications. Default = 1", type=int, default=1)
	parser_args.add_argument("-dumax", "--duplicationmaxlength", help="Maximum length of duplications. Default = 2", type=int, default=2)
	parser_args.add_argument("-dub", "--duplicationblock", help="Amount of bases blocked after duplication. Default = 1", type=int, default=1)
	parser_args.add_argument("-tl", "--translocation", help="Translocation rate. Default = 0", type=float, default=0)
	parser_args.add_argument("-tlmin", "--translocationminlength", help="Minimum length of translocations. Default = 1", type=int, default=1)
	parser_args.add_argument("-tlmax", "--translocationmaxlength", help="Maximum length of translocations. Default = 2", type=int, default=2)
	parser_args.add_argument("-tlb", "--translocationblock", help="Amount of bases blocked after translocations. Default = 1", type=int, default=1)
	parser_args.add_argument("-a", "--assembly", help="Assembly name for the VCF file. Default = 'Unknown'", default="Unknown")
	parser_args.add_argument("-s", "--species", help="Species name for the VCF file. Default = 'Unknown'", default="Unknown")
	parser_args.add_argument("-n", "--sample", help="Sample name for the VCF file", default="SAMPLE")
	parser_it = subparsers.add_parser("it", help="Generate interchromosomal translocations via commandline")
	parser_it.add_argument("interchromosomalrate", help="Rate of interchromosomal translocations.", type=float)
	args = parser.parse_args()
	filename = args.file
	outfile_basename = args.output
	rmt_file = None
	it_rate = None
	mut_rates = {}
	mut_lengs = {}
	mut_block={}
	titv = 1
	ignore_warnings=False
	assembly_name="Unknown"
	species_name="Unknown"
	sample_name="SAMPLE"
	if not hasattr(args, "rmtfile") and not hasattr(args, "interchromosomalrate"):
		mut_rates = {"sn": args.snp, "in": args.insert, "de": args.deletion, "iv": args.inversion,
					 "du": args.duplication, "tl": args.translocation}
		mut_lengs = {"min": {"in": args.insertminlength, "de": args.deletionminlength, "iv": args.inversionminlength, "du": args.duplicationminlength, "tl": args.translocationminlength}, "max":{"in": args.insertmaxlength, "de": args.deletionmaxlength, "iv": args.inversionmaxlength, "du": args.duplicationmaxlength, "tl": args.translocationmaxlength}}
		for typ in ["in", "de", "iv", "du", "tl"]:
			if mut_lengs["min"][typ] > mut_lengs["max"][typ]:
				print(f"ERROR: Minimum length for '{typ}' > maximum length.")
				return False, False, False, False, False, False, False, False, False, False, False, False
			if typ != "iv" and mut_lengs["min"][typ] < 1:
				print(f"ERROR: Minimum length for '{typ}' < 1.")
				return False, False, False, False, False, False, False, False, False, False, False, False
			if typ == "iv" and mut_lengs["min"][typ] < 2:
				print(f"ERROR: Minimum length for '{typ}' < 2.")
				return False, False, False, False, False, False, False, False, False, False, False, False
		mut_block = {"sn": args.snpblock, "in": args.insertblock, "de": args.deletionblock, "iv": args.inversionblock,
					 "du": args.duplicationblock, "tl": args.translocationblock}
		for key in mut_block.keys():
			if mut_block[key] < 1: mut_block[key] = 1 ; print(f"'{key}' block value was set to 1")
		assembly_name = args.assembly
		species_name = args.species
		sample_name = args.sample
		titv = args.transitionstransversions
		if not outfile_basename:
			outfile_basename = ".".join(filename.split(".")[:-1])+"_ms"
	elif hasattr(args, "interchromosomalrate"):
		it_rate = args.interchromosomalrate
		if not outfile_basename:
			outfile_basename = ".".join(filename.split(".")[:-1])+"_ms_it"
	else:
		rmt_file = args.rmtfile
		ignore_warnings=args.ignore_warnings
		if not outfile_basename:
			outfile_basename = ".".join(filename.split(".")[:-1])+"_ms_rmt"
	return filename, outfile_basename, rmt_file, mut_rates, mut_lengs, mut_block, assembly_name, species_name, sample_name, it_rate, ignore_warnings, titv


def save_mutations_vcf(vcf_name, ref_filename, fasta, chromosome, mut_list, assembly, species, sample, mode):
	"""
	Saves all mutations in a VCF file according to version 4.3.
	Parameters:
		vcf_name (str): Name or Path for the VCF file.
		ref_filename (str): Name of the reference Fasta.
		fasta (pyfaidx.Fasta): Chromosome information from any Fasta index file.
		chromosome (str): Name of the Chromosome.
		mut_list (list): A list containing sublists with all information about the inserted mutations.
		assembly (str): Name of the assembly.
		species (str): Name of the species.
		sample (str): Name of the sample.
		mode (str): Write / Append mode for the VCF writing.
	"""
	with open(vcf_name, mode) as hndl:
		if mode == "w":
			now = datetime.datetime.now()
			hndl.write("##fileformat=VCFv4.3\n")
			hndl.write("##filedate=" + str(now.year) + str(now.month) + str(now.day) + "\n")
			hndl.write("##source=Mutation-Simulator\n")
			hndl.write(f"##reference={ref_filename.split('/')[-1]}\n")
			for chrom in list(fasta.keys()):
				hndl.write("##contig=<ID=" + fasta[chrom].name + ",length=" + str(len(fasta[chrom])) + ",assembly=" + str(assembly) + ",species=\"" + str(species) + "\""">\n")
			hndl.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
			hndl.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
			hndl.write("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
			hndl.write("##ALT=<ID=INS,Description=\"Insert\">\n")
			hndl.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
			hndl.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
			hndl.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
			hndl.write("##ALT=<ID=DEL:ME,Description=\"Deletion of mobile element\">\n")
			hndl.write("##ALT=<ID=INS:ME,Description=\"Insertion of mobile element\">\n")
			hndl.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
			hndl.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n")
		too_long = []
		for i in tqdm(range(len(mut_list)), desc="Writing VCF", position=1, leave=False):
			entry = mut_list[i]
			#start, ref, alt, info = "", "", "", ""
			if entry[0] == "sn":
				start, ref, alt, info = entry[1] + 1, convert_ambiguous(fasta[chromosome][entry[1]]), convert_ambiguous(entry[2]), "."
			elif entry[0] == "in":
				if entry[1] > 0:
					start, ref, alt, info = entry[1], convert_ambiguous(fasta[chromosome][entry[1] - 1:entry[1]]), convert_ambiguous(fasta[chromosome][entry[1] - 1:entry[1]] + entry[3]), f"SVTYPE=INS;END={entry[1]};SVLEN={len(entry[3])}"
				else:
					start, ref, alt, info = entry[1] + 1, convert_ambiguous(fasta[chromosome][entry[1]]), convert_ambiguous(entry[3] + fasta[chromosome][entry[1]]), f"SVTYPE=INS;END={entry[1] + 1};SVLEN={len(entry[3])}"
			elif entry[0] == "du":
				start, ref, alt, info = entry[1] + 1, convert_ambiguous(fasta[chromosome][entry[1]:entry[2] + 1]), convert_ambiguous(fasta[chromosome][entry[1]:entry[2] + 1] + entry[3]), f"SVTYPE=DUP;END={entry[1]+len(fasta[chromosome][entry[1]:entry[2]+1])};SVLEN={len(entry[3])}"
			elif entry[0] == "de" or entry[0] == "tl":
				if entry[0] == "de":
					sv_type = "DEL"
				else:
					sv_type = "DEL:ME"
				if entry[1] > 0:
					start, ref, alt, info = entry[1], convert_ambiguous(fasta[chromosome][entry[1] - 1:entry[2] + 1]), convert_ambiguous(fasta[chromosome][entry[1] - 1:entry[1]]), f"SVTYPE={sv_type};END={entry[2] + 1};SVLEN=-{entry[2] - entry[1] + 1}"
				else:
					start, ref, alt, info = entry[1] + 1, convert_ambiguous(fasta[chromosome][0:entry[2] + 2]), convert_ambiguous(fasta[chromosome][entry[2] + 1]), f"SVTYPE={sv_type};END={entry[2] + 2};SVLEN=-{entry[2] - entry[1] + 1}"
			elif entry[0] == "iv":
				if str(fasta[chromosome][entry[1]:entry[2] + 1]) != str(entry[3][::-1]):
					start, ref, alt, info = entry[1] + 1, convert_ambiguous(fasta[chromosome][entry[1]:entry[2] + 1]), convert_ambiguous(entry[3][::-1]), f"SVTYPE=INV;END={entry[2] + 1};SVLEN=0"
				else:
					continue
			elif entry[0] == "tli":
				entry[1] = int(entry[1]) ###
				if not entry[4]:  # check for transloc inversion
					insert = str(entry[5])
				else:
					insert = str(entry[5][::-1])
				if entry[1] > 0:
					if entry[1] < len(fasta[chromosome]):
						start, ref, alt, info = entry[1], convert_ambiguous(fasta[chromosome][entry[1] - 1:entry[1]]), convert_ambiguous(fasta[chromosome][entry[1] - 1:entry[1]] + insert), f"SVTYPE=INS:ME;END={entry[1]};SVLEN={len(insert)}"
					else:
						too_long.append(entry)
						pbar.update(1)
						continue
				else:
					start, ref, alt, info = entry[1] + 1, convert_ambiguous(fasta[chromosome][entry[1]]), convert_ambiguous(insert + fasta[chromosome][entry[1]]), f"SVTYPE=INS:ME;END={entry[1] + 1};SVLEN={len(insert)}"

			hndl.write(f"{fasta[chromosome].name}\t{start}\t.\t{ref}\t{alt}\t.\t.\t{info}\tGT\t1\n")

		too_long = fix_too_long(too_long)  # creates a single vcf-entry for every entry past the last base
		if too_long:
			start, ref, alt, info = len(fasta[chromosome]), convert_ambiguous(fasta[chromosome][len(fasta[chromosome]) - 1]), convert_ambiguous(fasta[chromosome][len(fasta[chromosome]) - 1] + too_long), f"SVTYPE=INS:ME;END={len(fasta[chromosome]) + len(too_long)};SVLEN={len(too_long)}"
			hndl.write(f"{fasta[chromosome].name}\t{start}\t.\t{ref}\t{alt}\t.\t.\t{info}\tGT\t1\n")


def fix_too_long(tl):
	"""
	Takes a list of translocation entrys and returns them as a single entry.

	This is used in case of multiple translocation entries at the end of the chromosome.
	"""
	fixed = ""
	for entry in tl:
		if not entry[4]:
			fixed = f"{fixed}{entry[5]}"
		else:
			fixed = f"{fixed}{entry[5][::-1]}"
	return fixed


def convert_ambiguous(substring):
	"""Returns the inpur string with converted ambiguity codes"""
	return substring.replace("K", "G").replace("S", "C").replace("Y", "C").replace("M", "A").replace("W", "A")\
		.replace("R", "A").replace("B", "C").replace("D", "A").replace("H", "A").replace("V", "A").replace("-", "N")


def get_mutations(start, stop, mut_rates, mut_lengs, mut_block):
	"""
	Generates a list of mutations with translocations in a seperate list for a whole chromosome or RMT region.

	This function also returns a set of all blocked positions.

	Parameters:
		start (int): Startung position of the chromosome or region.
		stop (int): Stop position of the chromosome or region.
		mut_rates (dict): All rates for all mutation types with 2 letter acronymes.
		mut_lengs (dict): All mutation lengths settings for any mutation type.
		mut_block (dict): The range of the blocked area after any mutation for all types of mutations.

	Returns:
		mutations (list): A list containing all mutations except translocations.
		translocations (list): A list containing all translocational deletions.
		blocked_positions (set): Contains all positions which are blocked by mutations and mut_block.
	"""
	mut_type_chances, mut_rate = calc_mut_type_chances(mut_rates)
	if not mut_type_chances:
		return False, False, False
	mut_positions = get_mut_positions(start, stop, mut_rate)
	if not mut_positions:
		if type(mut_positions) == list:
			return -1, False, False
		else:
			return False, False, False
	mutations = []
	blocked_positions = blist([])
	translocations = []
	pbar = tqdm(total=len(mut_positions), desc=f"Finding mutations in range: {start+1}-{stop+1}", position=1, leave=False)
	i = 0
	mut_types = [choice(mut_type_chances[0], p=mut_type_chances[1], size=len(mut_positions))]
	while i < len(mut_positions):
		i_before=i
		mut = random_mutation_type(mut_positions[i], stop, mut_lengs, mut_types[0][i])
		if mut:
			mutations.append(mut)
			if not mut[0] in ["sn", "in"]:
				mut_range = range(mut[1], mut[2] + 1 + mut_block[mut[0]])
				if mut[0] == "tl":
					translocations.append(mut)
			else:
				mut_range = range(mut[1], mut[1] + 1 + mut_block[mut[0]])
			blocked_positions=blocked_positions+blist(mut_range)
			flag = True
			i += 1
			while flag and i < len(mut_positions):
				if not mut_positions[i] in mut_range:
					flag = False
				else:
					i += 1
		elif mut == None:
			i += 1
		elif mut == False:
			return False, False, False
		pbar.update(i-i_before)
	pbar.close()
	blocked_positions=array(blocked_positions)
	return mutations, translocations, blocked_positions


def calc_mut_type_chances(mut_rates):
	"""
	Converts a mut_rates dictionary to a single rate and returns the chances for every mutation type.

	Parameters:
		mut_rates (dict): All rates for all mutation types with 2 letter acronblist(fasta[:])ymes.

	Returns:
		mut_type_chances (dict): Contains all chances for any mutation type according to rate_sum.
		rate_sum (float): The combined mutation rate of all mutation types.
	"""
	rate_sum = 0
	for mut_rate in mut_rates:  # create sum for chances
		rate_sum = rate_sum + mut_rates[mut_rate]
	mut_type_chances = [[],[]]
	for mut_rate in mut_rates:  # calculate chance
		try:
			mut_type_chances[0].append(str(mut_rate))
			mut_type_chances[1].append(float(mut_rates[mut_rate] / rate_sum))
		except ZeroDivisionError:
			print("ERROR: Sum of mutation rates = 0")
			return False, False
	must_haves = ["sn", "in", "de", "du", "iv", "tl"]
	for mut_type in must_haves:  # add missing values
		if mut_type not in mut_type_chances[0]:
			mut_type_chances[0].append(mut_type)
			mut_type_chances[1].append(0.0)
	return mut_type_chances, rate_sum


def get_trans_inserts(translocations, data_length, blocked_positions, no_tl_regions):
	"""
	Creates an translocational insert for every translocational deletion.

	Parameters:
		translocations (list): A list containing all translocational deletions.
		data_length (int): Length of the chromosome.
		blocked_positions (set): Contains all positions which are blocked by mutations and mut_block.
		no_tl_regions (list): COntains all regions where translocations are blocked.

	Returns:
		trans_inserts (list): Contains all translocational inserts.
	"""
	if no_tl_regions:
		for region in no_tl_regions:
			blocked_positions.update(region)
	trans_inserts = []
	if len(translocations) <= (data_length - len(blocked_positions)):
		positions = choice(setdiff1d(arange(data_length), array(blocked_positions), True), len(translocations), replace=False)
		for i in tqdm(range(len(translocations)), desc="Finding transloc inserts", position = 1, leave=False):
			trans_inserts.append(["tli", positions[i], translocations[i][1], translocations[i][2],
								  transloc_invert((translocations[i][2] - translocations[i][1]) + 1)])
	else:
		positions = rnd.sample(set(range(0, data_length)) - set(blocked_positions), len(set(range(0, data_length)) - set(blocked_positions)))
		rnd.shuffle(translocations)
		for i in tqdm(range(len(translocations)), desc="Finding transloc inserts", position = 1, leave=False):
			if i < len(positions): # if i < len(set(range(0, data_length)) - set(blocked_positions)):
				trans_inserts.append(["tli", positions[i], translocations[i][1], translocations[i][2],
									  transloc_invert((translocations[i][2] - translocations[i][1]) + 1)])
			else:
				data_length = data_length + 1
				trans_inserts.append(["tli", data_length, translocations[i][1], translocations[i][2],
									  transloc_invert((translocations[i][2] - translocations[i][1]) + 1)])
	return trans_inserts


def get_mut_positions(start, stop, mut_rate):
	"""Returns a sorted set of all mutation starting points."""
	mut_count = int(((stop - start) + 1) * mut_rate)
	try:
		positions = sorted(choice(arange(start, stop + 1), mut_count, replace=False))
	except ValueError:
		print("ERROR: Rates too high.")
		positions = False
	return positions


def random_mutation_type(start, data_length, mut_lengs, mut_type):
	"""
	Assigns a stop value to a selected start value and mutation type.

	Parameters:
		start (int): Position of a mutation.
		data_length (int): Length of the chromosome or region.
		mut_lengs (dict): All mutation length settings for any mutation type.
		mut_type_chances (dict): Contains all chances for any mutation type with 2 letter acronymes.
		mut_type (str): 2 Letter acronym for the mutation type.

	Returns:
		mutation (list): A list containing a mutation type, a start value and if given a stop value.
	"""
	mutation = [mut_type]
	if mutation[0] not in mut_lengs["min"].keys() or mutation[0] not in mut_lengs["max"].keys():
		if not mutation[0]=="sn":
			print(f"ERROR: {mutation[0]} min/max length not defined")
			return False
	mutation.append(int(start))  # appends the start position to the mutation type list
	if mutation[0] == "iv":  # needs a minimum length of 2
		if start + mut_lengs["max"]["iv"] >= data_length:
			return None  # if it cant fit at the end it will be skipped
		stop = rnd.randint(start + mut_lengs["min"]["iv"] - 1, start + (mut_lengs["max"]["iv"] - 1))
		mutation.append(stop)
	elif mutation[0] == "in":
		stop = rnd.randint(start + mut_lengs["min"][mutation[0]] - 1, start + mut_lengs["max"][mutation[0]] - 1)
		mutation.append(stop)
	elif mutation[0] == "du":
		stop = rnd.randint(start + mut_lengs["min"][mutation[0]] - 1, start + mut_lengs["max"][mutation[0]] - 1)
		if stop > data_length:
			stop = data_length
		mutation.append(stop)  # appends the stop position
	elif mutation[0] in ["tl", "de"]:
		stop = rnd.randint(start + mut_lengs["min"][mutation[0]] - 1, start + mut_lengs["max"][mutation[0]] - 1)
		if stop > data_length:
			stop = data_length
		mutation.append(stop)  # appends the stop position
	return mutation


def transloc_invert(transloc_len):
	"""
	Picks True or False randomly if the length of the sequence is >1.

	This is used to determine if a translocation will invert or not.
	"""
	if rnd.randint(0, 1) == 0 or transloc_len <= 1:
		return False
	else:
		return True


def split_chrom(record, positions):
	"""
	Splits a chromosome in blist format at given positions and returns the splits in a list of blists.
	Parameters:
		record (blist): Chromosome.
		positions (list): List of the split-positions.
	Returns:
		splitted (list): List of blists containing the splitted chromosome.
	"""
	splitted=[]
	for i in range(len(positions)):
		if i ==0:
			splitted.append(record[:positions[i]])
		else:
			splitted.append(record[positions[i-1]:positions[i]])
		if i == len(positions)-1:
			splitted.append(record[positions[i]:])
	return splitted


def interchromosomal_transloc(fasta, rate, rmt):
	"""
	Creates interchromosomal translocations using breakends.
	Parameters:
		fasta (pyfaidx.Fasta): Chromosome information from any Fasta index file.
		rate (float): The rate of interchromosomal translocations.
		rmt (dict): The RMT information from load_rmt(). Can be None.
	Returns:
		records (list): Chromosome information in its mutated state.
		changes (dict): Contains lists with all breakends for every chromosome.
	"""
	blocked_choromosomes = []
	changes={}
	records=[]
	partners={}
	for key in fasta.keys():
		records.append(blist(fasta[key][:]))
	for i in range(len(records)): changes[i] = []
	if rmt:
		for i in range(len(rmt)):
			if rmt[i][1] == None:
				blocked_choromosomes.append(i)
	avail_chr=list(range(len(records)))
	for blocked_chr in blocked_choromosomes: # remove blocked from available
		avail_chr.remove(blocked_chr)
	for chrom in avail_chr: # find partner chromosome
		avail_chr.remove(chrom)
		if len(avail_chr)>0:
			partner = rnd.choice(avail_chr)
			partners[chrom]=partner
			avail_chr.remove(partner)
	chr_to_cross = list(partners.keys())
	if len(chr_to_cross)==0:
		print("ERROR: Not enough unblocked sequencess to simulate interchromosomal translocations. You need at least 2 sequences to create an interchromosomal translocation.")
		return False, False
	for chrom in chr_to_cross: #apply interchromosomal translocations
		sys.stdout.write(f"Generating interchromosomal translocations for sequence {chrom+1} and {partners[chrom]+1}\n")
		if rmt:
			amount = int(((len(records[chrom])+len(records[partners[chrom]])-4)/2) * ((rmt[chrom][-1]+rmt[partners[chrom]][-1])/2))
		else:
			amount = int((len(records[chrom])+len(records[partners[chrom]])-4)/2 * rate)
		try:
			pos_chrom = sorted(sample_with_minimum_distance(len(records[chrom]), amount, 2)) # keep at least 1 distance, cannot be 0 or last index
			pos_partner = sorted(sample_with_minimum_distance(len(records[partners[chrom]]), amount, 2))
		except ValueError:
			print("ERROR: Interchromosomal translocation rate too high.")
			return False, False
		splitted_chrom = split_chrom(records[chrom], pos_chrom)
		splitted_partner = split_chrom(records[partners[chrom]], pos_partner)
		for i in range(1, len(splitted_chrom), 2):# switching splits
			temp = splitted_chrom[i]
			splitted_chrom[i] = splitted_partner[i]
			splitted_partner[i] = temp
		records[chrom] = list(itertools.chain.from_iterable(splitted_chrom))
		records[partners[chrom]] = list(itertools.chain.from_iterable(splitted_partner))
		changes[chrom]={"self": pos_chrom, "partner": {"name": partners[chrom],"pos": pos_partner}}
		changes[partner]={"self": pos_partner, "partner": {"name": chrom,"pos": pos_chrom}}
	return records, changes


def save_it_bedpe(filename, fasta, it_changes):
	"""
	Saves all interchromosomal translocations in a BEDPE file.

	Parameters:
		filename (str): Name for the BEDPE file.
		fasta (pyfaidx.Fasta): Chromosome information from any Fasta index file.
		it_changes (dict): Dictionary containing all interchromosomal translocations for every chromosome.
	"""
	chromosome_names = list(fasta.keys())
	with open(filename, "w") as hndl:
		hndl.write("#chrom1\tstart1\tstop1\tchrom2\tstart2\tstop2\n")
		for chrom in list(it_changes.keys()):
			for i in range(0,len(it_changes[chrom]["self"]),2):
				if i != len(it_changes[chrom]["self"])-1:
					hndl.write(f'{chromosome_names[chrom]}\t{it_changes[chrom]["self"][i]}\t{it_changes[chrom]["self"][i+1]}\t{chromosome_names[it_changes[chrom]["partner"]["name"]]}\t{it_changes[chrom]["partner"]["pos"][i]}\t{it_changes[chrom]["partner"]["pos"][i+1]}\n')
				else:
					if not len(it_changes[chrom]["self"])%2 == 0:
						hndl.write(f'{chromosome_names[chrom]}\t{it_changes[chrom]["self"][i]}\t{len(fasta[chrom])}\t{chromosome_names[it_changes[chrom]["partner"]["name"]]}\t{it_changes[chrom]["partner"]["pos"][i]}\t{len(fasta[it_changes[chrom]["partner"]["name"]])}\n')
	return


def sublists_empty(lists):
	"""Returns True if every sublist of a list is empty."""
	return all(map(sublists_empty, lists)) if isinstance(lists, blist) else False


def subdicts_empty(dict):
	"""Returns True if every 'self' dict key in subdicts is empty."""
	is_empty = True
	for subdict in dict:
		if dict[subdict]["self"]: is_empty = False
	return is_empty


def load_rmt(rmt_file, fai, filename, ignore_warnings):
	"""
	Loads the given rmt file.

	Parameters:
		rmt_file (str): Name of the rmt file.
		fai (pyfaidx.Faidx): Fasta index file content.
		filename (str): Name of the fasta file.
		ignore_warnings (bool): Ignores RMT warnings.

	Returns:
		range_definitions (list): Contains the important rmt information.
		mut_flag (bool): Reflects if mutation information is in the range_definitions.
		it_flag (bool): Reflects if interchromosomal translocation information is in the range_definitions.
		species_name (str): Name of the species.
		assembly_name (str): Name of the assembly.
		sample_name (str): Name of the sample.
		mut_block (dict): The range of the blocked area after any mutation for all types of mutations.
		titv (float): Transition / Transversion ratio.
	"""
	print("Loading RMT")
	rmt = read_rmt(rmt_file)
	if rmt:
		range_definitions, rates_std, lengs_std, it_std, mut_flag, it_flag, meta= parse_rmt(rmt)
		passed, species_name, assembly_name, sample_name, mut_block, titv=rmt_meta_check(meta, filename, ignore_warnings)
		if not check_lengs_std(rates_std, lengs_std):
			passed = False
		if not check_lengs_rd(range_definitions):
			passed = False
		if not check_chr_exists(range_definitions, fai.index):
			passed = False
		if passed:
			range_definitions = set_missing_chr_2_std(range_definitions, len(fai.index), rates_std, lengs_std)
			range_definitions = add_missing_its(range_definitions, it_std)
			range_definitions = convert_ranges(range_definitions, fai)
			range_definitions = add_missing_ranges_as_std(range_definitions, fai, rates_std, lengs_std)
		else:
			return None, None, None, None, None, None, None, None
	else:
		return None, None, None, None, None, None, None, None
	return range_definitions, mut_flag, it_flag, species_name, assembly_name, sample_name, mut_block, titv


def check_lengs_std(rates_std, lengs_std):
	for typ in rates_std[0].keys():
		if typ != "sn":
			if typ not in lengs_std["min"]:
				print(f"ERROR: Minimum length for '{typ}' is missing in std.")
				return False
			if typ not in lengs_std["max"]:
				print(f"ERROR: Maximum length for '{typ}' is missing in std.")
				return False
			if lengs_std["min"][typ] > lengs_std["max"][typ]:
				print(f"ERROR: Minimum length for '{typ}' > maximum length in std.")
				return False
			if typ != "iv" and lengs_std["min"][typ] < 1:
				print(f"ERROR: Minimum length for '{typ}' < 1 in std.")
				return False
			if typ == "iv" and lengs_std["min"][typ] < 2:
				print(f"ERROR: Minimum length for '{typ}' < 2 in std.")
				return False
	return True


def check_lengs_rd(rd):
	for chr in rd:
		for rng in chr[1]:
			if rng[1]:
				for typ in rng[1].keys():
					if typ !="sn":
						if typ not in rng[2]["min"]:
							print(f"ERROR: Minimum length for '{typ}' is missing in range {rng[0]} on chr {chr[0]}.")
							return False
						if typ not in rng[2]["max"]:
							print(f"ERROR: Maximum length for '{typ}' is missing in range {rng[0]} on chr {chr[0]}.")
							return False
						if rng[2]["min"][typ] > rng[2]["max"][typ]:
							print(f"ERROR: Minimum length for '{typ}' > maximum length in range {rng[0]} on chr {chr[0]}.")
							return False
						if typ != "iv" and rng[2]["min"][typ] < 1:
							print(f"ERROR: Minimum length for '{typ}' < 1 in range {rng[0]} on chr {chr[0]}.")
							return False
						if typ == "iv" and rng[2]["min"][typ] < 2:
							print(f"ERROR: Minimum length for '{typ}' < 2 in range {rng[0]} on chr {chr[0]}.")
							return False
	return True


def check_chr_exists(rd, fai):
	chroms = list(fai.keys())
	for chr in rd:
		if len(chroms) < int(chr[0]):
			print(f"ERROR: Fasta only has {len(chroms)} sequence/s. Chr index {chr[0]} in RMT is out of range.")
			return False
	return True


def read_rmt(file):
	"""Reads an RMT file."""
	list_of_lines_raw = []
	list_of_lines = []
	try:
		with open(file, "r") as hndl:
			for line in hndl.readlines():
				if not line.startswith("#"):
					list_of_lines_raw.append(line.strip().split("#")[0])
			for line in list_of_lines_raw:
				if not line == "":
					list_of_lines.append(line)
	except FileNotFoundError:
		print("ERROR: Rmt file can not be found")
		list_of_lines = None
	return list_of_lines


def parse_rmt(rmt):
	"""Parses an RMT file."""
	meta={}
	rates_std = [{}, False]
	lengs_std = {"min": {}, "max": {}}
	range_definitions = []
	it_std = False
	rate_indicators = ["sn", "in", "de", "iv", "du", "tl"]
	max_leng_indicators = ["inmax", "demax", "ivmax", "dumax", "tlmax"]
	min_leng_indicators = ["inmin", "demin", "ivmin", "dumin", "tlmin"]
	flag = 0
	mut_flag=False
	it_flag=False
	for i in range(0, len(rmt)):
		line = rmt[i].split(" ")
		if flag == 0 and not line[0].lower() == "std":
			for info in line:
				info = info.split("=")
				try:
					meta[info[0]] = " ".join(info[1:]+line[1:])
				except IndexError:
					pass
		elif flag == 1 and not line[0].lower() == "chr":
			for n in range(0, len(line)):
				if line[n].lower() in rate_indicators:
					mut_flag=True
					if line[n].lower() == "tl":
						rates_std[1] = True
						rates_std[0][line[n]] = float(line[n + 1]) / 2  # half the tl rate
					else:
						rates_std[0][line[n]] = float(line[n + 1])
				elif line[n].lower() in max_leng_indicators:
					lengs_std["max"][line[n][:2]] = int(line[n + 1])
				elif line[n].lower() in min_leng_indicators:
					lengs_std["min"][line[n][:2]] = int(line[n + 1])
				elif line[n].lower() == "it":
					if not line[n + 1].lower() == "none":
						it_flag=True
						it_std = float(line[n + 1]) # it rate gets normalized in interchromosomal_transloc()
					else:
						it_std = None
				elif line[n].lower() == "None":
					rates_std = None
					lengs_std = None
		elif flag == 2 and not line[0].lower() == "chr":
			if not line[0].startswith("it"):
				range_definitions[-1][1].append([line[0], {}, {"min": {}, "max": {}}, False])  # chr mutrates mutlen tlinsertsflag
				for n in range(0, len(line)):
					if line[n].lower() in rate_indicators:
						mut_flag=True
						if line[n].lower() == "tl":
							range_definitions[-1][1][-1][3] = True  # sets a flag for containing tl
							range_definitions[-1][1][-1][1][line[n]] = float(line[n + 1]) / 2
						else:
							range_definitions[-1][1][-1][1][line[n]] = float(line[n + 1])
					elif line[n].lower() in min_leng_indicators:
						range_definitions[-1][1][-1][2]["min"][line[n][:2]] = float(line[n + 1])
					elif line[n].lower() in max_leng_indicators:
						range_definitions[-1][1][-1][2]["max"][line[n][:2]] = float(line[n + 1])
					elif line[n].lower() == "none":
						range_definitions[-1][1][-1][1] = None
						range_definitions[-1][1][-1][2] = None
						range_definitions[-1][1][-1][3] = None
			elif line[0].lower() == "it":
				if not line[1].lower() == "none":
					it_flag = True
					range_definitions[-1].append(float(line[1]))
				else:
					range_definitions[-1].append(None)
		if line[0].lower() == "std":
			flag = 1
		if line[0].lower() == "chr":
			range_definitions.append([line[1], []])
			flag = 2
	return range_definitions, rates_std, lengs_std, it_std, mut_flag, it_flag, meta


def rmt_meta_check(meta,filename,ignore_warnings):
	"""Checks the RMT meta information."""
	check=True
	mut_indicator=["sn", "in", "de", "du", "iv", "tl"]
	mut_block={}
	if not ignore_warnings:
		if "fasta" in meta.keys():
			if not meta["fasta"]==filename.split("/")[-1]:
				yn = input("Fastaname does not match rmt, ignore? [y/n]")
				if not yn.lower() == "y": check = False
		if "md5" in meta.keys():
			if not get_md5(filename)==meta["md5"]:
				yn = input("The fasta's md5-hash does not match rmt, ignore? [y/n]")
				if not yn.lower() == "y":
					check = False
	if "species_name" in meta.keys():
		species_name = meta["species_name"]
	else:
		species_name = "Unknown"
	if "assembly_name" in meta.keys():
		assembly_name = meta["assembly_name"]
	else:
		assembly_name = "Unknown"
	if "sample_name" in meta.keys():
		sample_name = meta["sample_name"]
	else:
		sample_name = "SAMPLE"
	if "titv" in meta.keys():
		titv = meta["titv"]
	else:
		titv = 1
	for key in meta.keys():
		for indicator in mut_indicator:
			if key == indicator+"_block":
				mut_block[indicator] = int(meta[key])
	for indicator in mut_indicator:
		if indicator not in mut_block.keys():
			mut_block[indicator] = 0
	return check, species_name, assembly_name, sample_name, mut_block, titv


def get_md5(filename):
	"""Calculated the md5 hash of a file."""
	hash_md5 = hashlib.md5()
	with open(filename, "rb") as f:
		for chunk in iter(lambda: f.read(4096), b""):
			hash_md5.update(chunk)
	return hash_md5.hexdigest()


def set_missing_chr_2_std(rd, chr_count, rates_std, lengs_std):
	"""Sets all mssing chromosomes to standard in an RMT."""
	chr_contigs = []
	for i in range(0, len(rd)):
		chr_contigs.append(int(rd[i][0]) - 1)
		rd[i].pop(0)
	chr_contigs = set(range(0, chr_count)) - set(chr_contigs)  # missing chrs
	for chromosome in chr_contigs:
		rd.insert(chromosome, [[["1-END",rates_std[0], lengs_std, rates_std[1]]]])
	return rd


def add_missing_its(rd, it_std):
	"""Adds missing interchromosomal translocation rates in the rd section."""
	for i in range(len(rd)):
		if len(rd[i]) < 2:
			rd[i].append(it_std)
	return rd


def convert_ranges(rd, fai):
	"""Converts all rd section ranges to python readable."""
	keys=list(fai.index.keys())
	for i in range(0, len(rd)):
		for n in range(len(rd[i][0])):
			converted_range = rd[i][0][n][0].split("-")
			converted_range[0] = int(converted_range[0]) - 1
			if converted_range[1] == "END":
				converted_range[1] = fai.index[keys[i]].rlen
			converted_range[1] = int(converted_range[1]) - 1
			rd[i][0][n][0] = [converted_range[0], converted_range[1]]
	return rd


def add_missing_ranges_as_std(rd, fai, rates_std, lengs_std):
	"""Adds all missing ranges and sets them to standard in an RMT."""
	keys = list(fai.index.keys())
	for i in range(0, len(rd)):
		if not rd[i][0]:
			rd[i][0].append([[0, fai.index[keys[i]].rlen - 1], rates_std[0], lengs_std, rates_std[1]])
		else:
			if not rd[i][0][0][0][0] == 0:
				rd[i][0].insert(0, [[0, rd[i][0][0][0][0] - 1], rates_std[0], lengs_std, rates_std[1]])
			n = 0
			while n < len(rd[i][0]):
				if n < len(rd[i][0]) - 1:
					if not rd[i][0][n][0][1] + 1 == rd[i][0][n + 1][0][0]:
						rd[i][0].insert(n + 1, [[rd[i][0][n][0][1] + 1, rd[i][0][n + 1][0][0] - 1], rates_std[0], lengs_std, rates_std[1]])
						n += 1
				else:
					if rd[i][0][n][0][1] < fai.index[keys[i]].rlen - 1:
						rd[i][0].insert(n + 1, [[rd[i][0][n][0][1] + 1, fai.index[keys[i]].rlen - 1], rates_std[0], lengs_std,
											 rates_std[1]])
						n += 1
				n += 1
	return rd


def mutate(fasta, mut_list, titv):
	"""
	This function creates the simulated sequence from a reference sequence and a list of mutations.

	Parameters:
		fasta (pyfaidx.Fasta): Chromosome information from any Fasta index file.
		mut_list (list): A list containing sublists with all information about the inserted mutations.
		titv (float): Transition / Transversion ratio.
	"""
	data = blist(fasta[:])
	for i in tqdm(range(len(mut_list)), desc="Mutating", position=1, leave=False):
		if mut_list[i][0] == "sn":
			alt = get_snp(data[mut_list[i][1]], titv)
			data[mut_list[i][1]]=alt
			mut_list[i].append(alt)  # appends the altered base, this is needed for the vcf file
		elif mut_list[i][0] == "in":
			insert = get_insert((mut_list[i][2]+1 -mut_list[i][1]))
			data=data[:mut_list[i][1]]+blist(insert)+data[mut_list[i][1]:]
			mut_list[i].append(insert)
		elif mut_list[i][0] == "de" or mut_list[i][0] == "tl":
			del(data[mut_list[i][1]:mut_list[i][2]+1])
		elif mut_list[i][0] == "iv":
			invert = blist(data[mut_list[i][1]:mut_list[i][2] + 1])
			del(data[mut_list[i][1]:mut_list[i][2] + 1])
			data=data[:mut_list[i][1]] + invert[::-1] +data[mut_list[i][1]:]
			mut_list[i].append("".join(invert))
		elif mut_list[i][0] == "du":
			dupe = blist(data[mut_list[i][1]:mut_list[i][2] + 1])
			data=data[:mut_list[i][2] + 1] + dupe + data[mut_list[i][2] + 1:]
			mut_list[i].append("".join(dupe))
		elif mut_list[i][0] == "tli":
			insert = blist(fasta[mut_list[i][2]:mut_list[i][3] + 1])
			if mut_list[i][4]:
				data = data[:mut_list[i][1]] + insert[::-1] + data[mut_list[i][1]:]
			else:
				data = data[:mut_list[i][1]] + insert + data[mut_list[i][1]:]
			mut_list[i].append("".join(insert))
	return data, mut_list


def get_snp(base, ti_tv):
	"""Returns a mutated base with the probability given by the transition to transversion rate ti_tv."""
	pTi = ti_tv * (1 / (ti_tv + 1))
	p = random()
	if p <= pTi:
		return get_ti_Base(base)
	else:
		return get_tv_Base(base)


def get_tv_Base(base):
	"""Returns one of the two transversion bases with equal probability."""
	if base not in ["A", "T", "C", "G"]: base = convert_ambiguous(base)
	return {"A": ["T", "C"], "G": ["C", "T"], "T": ["G", "A"], "C": ["A", "G"], "N": ["N", "N"]}[base][randint(2)] # replace randint


def get_ti_Base(base):
	"""Returns the transition base."""
	if base not in ["A", "T", "C", "G"]: base = convert_ambiguous(base)
	return {"A": "G", "G": "A", "T": "C", "C": "T", "N": "N"}[base]


def get_insert(leng):
	"""Returns a random sequence of specified length."""
	return "".join(choice(["A", "T", "G", "C"], leng))


def write_fasta(filename, chromosome, header, linebases, mode):
	"""
	Saves a list of strings in a Fasta file.

	Parameters:
		filename (str): Name of the output Fasta file.
		chromosome (list): Chromosome sequence information.
		header (str): Fasta sequence header.
		linebases (int): Number of bases per line.
		mode (str): Write / Append mode.
	"""
	with open(filename, mode) as hndl:
		hndl.write(f">{header}\n")
		for i in tqdm(range(0,len(chromosome),linebases), desc="Writing Sequence", position=1, leave=False):
			hndl.write("".join(chromosome[i:i+linebases])+"\n")
	return


def sample_with_minimum_distance(n=40, k=4, d=10):
	"""
	Sample of k elements from range(n), with a minimum distance d. Cannot be 0 or n.
	"""
	sample = rnd.sample(range(1, n-(k-1)*(d-1)), k)
	indices = sorted(range(len(sample)), key=lambda i: sample[i])
	return [s + (d-1)*r for s, r in zip(sample, sorted(indices, key=lambda i: indices[i]))]


if __name__ == "__main__":
	try:
		main()
	except KeyboardInterrupt:
		print("\n KeyboardInterrupt")
		sys.exit(1)
