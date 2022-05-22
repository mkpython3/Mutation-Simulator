from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path
from typing import TYPE_CHECKING

from ._version import __version__
from .defaults import Defaults

if TYPE_CHECKING:
	from argparse import Namespace


def add_outfile_names(args: Namespace) -> Namespace:
	"""Adds oufasta, outfastait, outvcf and outbedpe based on outbase.
	:param args: Commandline arguments
	:return args: Commandline arguments with additions
	:note: Outbase will be modified
	"""
	try:
		args.outbase = args.outbase.with_stem(args.outbase.stem + "_ms")
	except ValueError:
		args.outbase = args.outbase / (args.infile.stem + "_ms")
	args.outfasta = args.outbase.with_suffix(args.infile.suffix)
	args.outfastait = args.outfasta.with_stem(args.outfasta.stem + "_it")
	args.outvcf = args.outfasta.with_suffix(".vcf")
	args.outbedpe = args.outfastait.with_suffix(".bedpe")
	return args


def get_args() -> "Namespace":
	"""Returns commandline arguments.
	:return args: Commandline arguments
	"""
	parser = ArgumentParser(prog="mutation-simulator",
			description=
			"See https://github.com/mkpython3/Mutation-Simulator for more information about this program."
							)
	parser.add_argument("infile",
			type=Path,
			help="Path of the reference Fasta file")
	parser.add_argument("-o",
			"--output",
			type=Path,
			help="Path/Basename for the output files (without file extension)",
			default=Defaults.OUTBASE,
			dest="outbase")
	parser.add_argument("--ignore-warnings",
			help="Silences warnings",
			action='store_true',
			default=Defaults.IGNORE_WARNINGS)
	parser.add_argument("-v",
			"--version",
			action="version",
			version=f"Mutation-Simulator {__version__}")

	subparsers = parser.add_subparsers(dest="mode",
			help=
			"Generate mutations or interchromosomal translocations via RMT or arguments"
										)
	subparsers.required = True

	parser_rmt = subparsers.add_parser("rmt",
			help="Use random mutation table instead of arguments")
	parser_rmt.add_argument("rmtfile", type=Path, help="Path to the RMT file")

	parser_args = subparsers.add_parser("args",
			help="Use commandline arguments for mutations instead of RMT")
	parser_args.add_argument("-sn",
			"--snp",
			help=f"SNP rate. Default = {Defaults.RATE}",
			type=float,
			default=Defaults.RATE)
	parser_args.add_argument("-snb",
			"--snpblock",
			help=
			f"Amount of bases blocked after SNP. Default = {Defaults.BLOCK}",
			type=int,
			default=Defaults.BLOCK)
	parser_args.add_argument("-titv",
			"--transitionstransversions",
			help=
			f"Ratio of transitions:transversions likelihood. Default = {Defaults.TITV}",
			type=float,
			default=Defaults.TITV)
	parser_args.add_argument("-in",
			"--insert",
			help=f"Insert rate. Default = {Defaults.RATE}",
			type=float,
			default=Defaults.RATE)
	parser_args.add_argument("-inmin",
			"--insertminlength",
			help=f"Minimum length of inserts. Default = {Defaults.MINLEN}",
			type=int,
			default=Defaults.MINLEN)
	parser_args.add_argument("-inmax",
			"--insertmaxlength",
			help=f"Maximum length of inserts. Default = {Defaults.MAXLEN}",
			type=int,
			default=Defaults.MAXLEN)
	parser_args.add_argument("-inb",
			"--insertblock",
			help=
			f"Amount of bases blocked after insert. Default = {Defaults.BLOCK}",
			type=int,
			default=Defaults.BLOCK)
	parser_args.add_argument("-de",
			"--deletion",
			help=f"Deletion rate. Default = {Defaults.RATE}",
			type=float,
			default=Defaults.RATE)
	parser_args.add_argument("-demin",
			"--deletionminlength",
			help=f"Minimum length of deletions. Default = {Defaults.MINLEN}",
			type=int,
			default=Defaults.MINLEN)
	parser_args.add_argument("-demax",
			"--deletionmaxlength",
			help=f"Maximum length of deletions. Default = {Defaults.MAXLEN}",
			type=int,
			default=Defaults.MAXLEN)
	parser_args.add_argument("-deb",
			"--deletionblock",
			help=
			f"Amount of bases blocked after deletion. Default = {Defaults.BLOCK}",
			type=int,
			default=Defaults.BLOCK)
	parser_args.add_argument("-iv",
			"--inversion",
			help=f"Inversion rate. Default = {Defaults.RATE}",
			type=float,
			default=Defaults.RATE)
	parser_args.add_argument("-ivmin",
			"--inversionminlength",
			help=f"Minimum length of inversion. Default = {Defaults.IV_MINLEN}",
			type=int,
			default=Defaults.IV_MINLEN)
	parser_args.add_argument("-ivmax",
			"--inversionmaxlength",
			help=f"Maximum length of inversion. Default = {Defaults.IV_MAXLEN}",
			type=int,
			default=Defaults.IV_MAXLEN)
	parser_args.add_argument("-ivb",
			"--inversionblock",
			help=
			f"Amount of bases blocked after inversion. Default = {Defaults.BLOCK}",
			type=int,
			default=Defaults.BLOCK)
	parser_args.add_argument("-du",
			"--duplication",
			help=f"Duplication rate. Default = {Defaults.RATE}",
			type=float,
			default=Defaults.RATE)
	parser_args.add_argument("-dumin",
			"--duplicationminlength",
			help=f"Minimum length of duplications. Default = {Defaults.MINLEN}",
			type=int,
			default=Defaults.MINLEN)
	parser_args.add_argument("-dumax",
			"--duplicationmaxlength",
			help=f"Maximum length of duplications. Default = {Defaults.MAXLEN}",
			type=int,
			default=Defaults.MAXLEN)
	parser_args.add_argument("-dub",
			"--duplicationblock",
			help=
			f"Amount of bases blocked after duplication. Default = {Defaults.BLOCK}",
			type=int,
			default=Defaults.BLOCK)
	parser_args.add_argument("-tl",
			"--translocation",
			help=f"Translocation rate. Default = {Defaults.RATE}",
			type=float,
			default=Defaults.RATE)
	parser_args.add_argument("-tlmin",
			"--translocationminlength",
			help=
			f"Minimum length of translocations. Default = {Defaults.MINLEN}",
			type=int,
			default=Defaults.MINLEN)
	parser_args.add_argument("-tlmax",
			"--translocationmaxlength",
			help=
			f"Maximum length of translocations. Default = {Defaults.MAXLEN}",
			type=int,
			default=Defaults.MAXLEN)
	parser_args.add_argument("-tlb",
			"--translocationblock",
			help=
			f"Amount of bases blocked after translocations. Default = {Defaults.BLOCK}",
			type=int,
			default=Defaults.BLOCK)
	parser_args.add_argument("-a",
			"--assembly",
			help=
			f"Assembly name for the VCF file. Default = '{Defaults.ASSEMBLY_NAME}'",
			default=Defaults.ASSEMBLY_NAME)
	parser_args.add_argument("-s",
			"--species",
			help=
			f"Species name for the VCF file. Default = '{Defaults.SPECIES_NAME}'",
			default=Defaults.SPECIES_NAME)
	parser_args.add_argument("-n",
			"--sample",
			help=
			f"Sample name for the VCF file. Default = '{Defaults.SAMPLE_NAME}'",
			default=Defaults.SAMPLE_NAME)

	parser_it = subparsers.add_parser("it",
			help="Generate interchromosomal translocations via the command line"
										)
	parser_it.add_argument("interchromosomalrate",
			help="Rate of interchromosomal translocations",
			type=float)
	args = parser.parse_args()

	return add_outfile_names(args)
