from __future__ import annotations

from datetime import datetime
from typing import TYPE_CHECKING

if TYPE_CHECKING:
	from pathlib import Path

	from pyfaidx import Fasta


class VcfWriterError(Exception):
	"""Raised when the writer can not write to a file."""


class VcfRecord:
	"""Contains VCF fields."""

	def __init__(self,
			svtype: str = "",
			start: int = 0,
			end: int = 0,
			len: int = 0,
			ref: str = "",
			alt: str = ""):
		"""Constructor.
		:param svtype: SV Type for the VCF
		:param start: Starting position
		:param end: End position
		:param len: SVLEN for the VCF
		:param ref: REF field
		:param alt: ALT field
		"""
		self.svtype = svtype
		self.start = start
		self.end = end
		self.len = len
		self.ref = ref
		self.alt = alt

	def __repr__(self) -> str:
		return f"{self.svtype} {self.start} {self.end} {self.len} {self.ref} {self.alt}"

	@property
	def info(self):
		"""Returns the info string for a VCF entry.
		:return: Info string
		"""
		if self.svtype != "sn":
			return f"SVTYPE={self.svtype};END={self.end};SVLEN={self.len}"
		else:
			return "."


class VcfWriter:
	"""Writes a new VCF file."""

	def __init__(self, fname: Path):
		"""Constructor.
		:param fname: VCF file to create
		"""
		try:
			self.__outf = open(fname, "w")
		except (IOError) as e:
			raise VcfWriterError(f"Cannot write to VCF file {fname} {e}")

	def __del__(self):
		self.close()

	def close(self):
		"""Closes the file handle."""
		self.__outf.close()

	def write_header(self, input_fasta: Path, fasta: Fasta, assembly_name: str,
			species_name: str, sample_name: str):
		"""Writes the VCF header.
		:param input_fasta: Name/Path of the reference Fasta
		:param fasta: Fasta file used in the simulation
		:param assembly_name: Name of the assembly
		:param species_name: Name of the species
		:param sample_name: Name for the sample column
		"""
		now = datetime.now()
		self.__outf.write("##fileformat=VCFv4.3\n")
		self.__outf.write(f"##filedate={now.year}{now.month}{now.day}\n")
		self.__outf.write("##source=Mutation-Simulator\n")
		self.__outf.write(f"##reference={input_fasta}\n")
		for chrom in list(fasta.keys()):
			self.__outf.write(
					f"##contig=<ID={fasta[chrom].name},length={len(fasta[chrom])},assembly={assembly_name},species=\"{species_name}\">\n"
			)
		self.__outf.write(
				"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
		)
		self.__outf.write(
				"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n"
		)
		self.__outf.write(
				"##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"
		)
		self.__outf.write("##ALT=<ID=INS,Description=\"Insert\">\n")
		self.__outf.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
		self.__outf.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
		self.__outf.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
		self.__outf.write(
				"##ALT=<ID=DEL:ME,Description=\"Deletion of mobile element\">\n"
		)
		self.__outf.write(
				"##ALT=<ID=INS:ME,Description=\"Insertion of mobile element\">\n"
		)
		self.__outf.write(
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
		)
		self.__outf.write(
				f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
		)

	def write(self, record: VcfRecord, seq_name: str):
		"""Writes a VCF entry.
		:param VcfRecord: Object to write as a VCF line
		:param seq_name: Sequence name for the chromosome column
		"""
		if record.ref != record.alt:
			self.__outf.write(
					f"{seq_name}\t{record.start}\t.\t{record.ref}\t{record.alt}\t.\t.\t{record.info}\tGT\t1\n"
			)
