from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
	from pathlib import Path


class FastaWriterError(Exception):
	"""Raised when the writer can not write to a file."""


class FastaWriter:
	"""Writes a new Fasta file."""

	def __init__(self, fname: Path):
		"""Constructor.
		:param fname: Fasta file to create
		"""
		try:
			self.__outf = open(fname, "w")
		except (IOError) as e:
			raise FastaWriterError(f"Cannot write to Fasta file {fname} {e}")
		self.__written = 0
		self.__bpl = 60

	def __del__(self):
		self.close()

	def close(self):
		"""Closes the file handle."""
		self.__outf.close()

	def set_bpl(self, bpl: int):
		"""Sets the bases per line limit for the Fasta file.
		:param bpl: Bases per line limit
		"""
		self.__bpl = bpl

	def write_header(self, header: str):
		"""Writes the fasta header.
		:param header: Fasta header (without ">")
		"""
		if self.__written != 0:
			self.__outf.write("\n")
		self.__outf.write(f">{header}\n")
		self.__written = 0

	def write(self, base: str):
		"""Writes a base to the Fasta.
		:param base: The base
		"""
		self.__outf.write(base)
		self.__written += 1

		if self.__written % self.__bpl == 0:
			self.__outf.write("\n")
			self.__written = 0

	def write_multi(self, bases: str | list[str]):
		"""Writes multiple base to the Fasta.
		:param bases: The bases
		"""
		for base in bases:
			self.write(base)
