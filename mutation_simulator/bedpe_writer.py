from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
	from pathlib import Path


class BedpeWriterError(Exception):
	"""Raised when the writer can not write to a file."""


class BedpeWriter:
	"""Writes a new BEDPE file."""

	def __init__(self, fname: Path):
		"""Constructor.
		:param fname: BEDPE file to create
		"""
		try:
			self.__outf = open(fname, "w")
		except (IOError) as e:
			raise BedpeWriterError(f"Cannot write to BEDPE file {fname} {e}")

	def __del__(self):
		self.close()

	def close(self):
		"""Closes the file handle."""
		self.__outf.close()

	def write_header(self):
		"""Writes the BEDPE header."""
		self.__outf.write("#chrom1\tstart1\tstop1\tchrom2\tstart2\tstop2\n")

	def write(self, chrom: str, bp_chrom: list[int], chrom_len_pre_it: int,
			partner: str, bp_partner: list[int], partner_len_pre_it: int):
		"""Writes multiple new BEDPE entries according to the breakpoints.
		:param chrom: Name of the first chromosome
		:param bp_chrom: List of breakpoints on the first chromosome
		:param chrom_len_pre_it: Length of the first chromosome before it
		:param partner: Name of the second chromosome
		:param bp_partner: List of breakpoints on the second chromosome
		:param partner_len_pre_it: Length of the second chromosome before it
		"""
		for i in range(0, len(bp_chrom), 2):
			if i != len(bp_chrom) - 1:
				self.__outf.write(
						f'{chrom}\t{bp_chrom[i]}\t{bp_chrom[i+1]}\t{partner}\t{bp_partner[i]}\t{bp_partner[i+1]}\n'
				)
			else:
				if not len(bp_chrom) % 2 == 0:
					self.__outf.write(
							f'{chrom}\t{bp_chrom[i]}\t{chrom_len_pre_it}\t{partner}\t{bp_partner[i]}\t{partner_len_pre_it}\n'
					)
