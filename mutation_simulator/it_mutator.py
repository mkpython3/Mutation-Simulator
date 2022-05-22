from __future__ import annotations

from random import choice, shuffle
from sys import stderr
from typing import TYPE_CHECKING, Tuple

from tqdm import trange

from .bedpe_writer import BedpeWriter
from .colors import Colors
from .fasta_writer import FastaWriter
from .util import pairwise, sample_with_minimum_distance

if TYPE_CHECKING:
	from argparse import Namespace

	from pyfaidx import Fasta

	from .rmt import SimulationSettings


class ITMutator:
	"""Stores required data for the simulation and performs interchromosomal
	translocations.
	"""

	def __init__(self, args: Namespace, fasta: Fasta, sim: SimulationSettings):
		"""Constructor.
		:param args: Commandline arguments
		:param fasta: Fasta file used in the simulation
		:param sim: Generated SimulationSettings from args, it or rmt mode
		"""
		self.__args = args
		self.__fasta = fasta
		self.__sim = sim

		self.__fasta_writer = FastaWriter(args.outfastait)
		self.__bedpe_writer = BedpeWriter(args.outbedpe)

		self.__find_avail_chroms()
		self.__assign_parters()

	def __del__(self):
		self.close()

	def close(self):
		"""Closes the fasta_writer and bedpe_writer filehandles. Does not close
		the fasta.
		"""
		self.__fasta_writer.close()
		self.__bedpe_writer.close()

	def __find_avail_chroms(self):
		"""Finds all chromosomes elligeble for interchromosomal translocations."""
		self.__avail_chroms = [
				chrom.number for chrom in self.__sim.chromosomes if
				chrom.it_rate is not None  # intentional is comparison, 0 is ok
				and len(self.__fasta[chrom.number]) > 2
		]

	def __assign_parters(self):
		"""Assigns a random unique partner to all chromosomes."""
		self.__partners = {}
		remain_chr = self.__avail_chroms
		shuffle(self.__avail_chroms)
		for chrom in self.__avail_chroms:
			remain_chr.remove(chrom)
			if remain_chr:
				partner = choice(remain_chr)
				self.__partners[partner] = chrom
				self.__partners[chrom] = partner
				remain_chr.remove(partner)

	def __get_pairs_to_mut_once(self) -> list[int]:
		"""Returns every chromosome pair to mutate only once. (1,2) will not
		get returned a second time as (2,1) which would normally be the case.
		:return chroms: List of chromosome indices
		"""
		chroms = list(self.__partners.keys())
		for chrom, partner in self.__partners.items():
			if chrom in chroms:
				chroms.remove(partner)
		return chroms

	def __get_bp_amount(self, seq_len1: int, rate1: float, seq_len2: int,
			rate2: float) -> int:
		"""Get the amount of breakpoints.
		:param seq_len1: Length of sequence 1
		:param rate1: It rate for sequence 1
		:param seq_len2: Length of sequence 2
		:param rate2: It rate for sequence 2
		:return breakpoint_amount: Number of breakpoints to generate
		"""
		return int((seq_len1 + seq_len2 - 4) / 2 * ((rate1 + rate2) / 2))

	def __get_breakpoints(self, chrom: int, seq_len1: int, seq_len2: int,
			bp_amount: int) -> Tuple[list[int], list[int]]:
		"""Generates the given amount of breakpoints on each sequence when
		possible.
		:param chrom: Chromosome index
		:param seq_len1: Length of sequence 1
		:param seq_len2: Length of sequence 2
		:param bp_amount: Number of breakpoints to generate
		:return bp_chrom: Breakpoints for sequence 1
		:return bp_partner: Breakpoints for sequence 2
		"""
		bp_chrom = []
		bp_partner = []

		try:
			# Keep 1 base at least inbetween breakpoints
			bp_chrom = sorted(
					sample_with_minimum_distance(seq_len1, bp_amount, 1))
			bp_partner = sorted(
					sample_with_minimum_distance(seq_len2, bp_amount, 1))
		except ValueError:
			if not self.__args.ignore_warnings:
				print(f"{Colors.warn}WARNING: Interchromosomal translocation rate too high for sequence {chrom+1} and {self.__partners[chrom]+1}.{Colors.norm}",
						file=stderr)
		return bp_chrom, bp_partner

	def __write_with_bp(self, chrom: str, bp_chrom: list[int], chrom_len: int,
			partner: str, bp_partner: list[int], partner_len: int):
		"""Writes a chromosome to the opened Fasta file with breakpoints.
		:param chrom: Name of the first chromosome
		:param bp_chrom: Breakpoints of the first chromosome
		:param chrom_len: Length of the first sequence
		:param partner: Name of the second chromosome
		:param bp_partner: Breakpoints of the second chromosome
		:param partner_len: Length of the second sequence
		"""
		bp_chrom = [0] + bp_chrom + [chrom_len]
		bp_partner = [0] + bp_partner + [partner_len]

		for i, (interval_chrom, interval_partner) in enumerate(
				zip(pairwise(bp_chrom), pairwise(bp_partner))):
			if i % 2:
				self.__fasta_writer.write_multi(
						self.__fasta.get_seq(partner, interval_partner[0] + 1,
						interval_partner[1]))
			else:
				self.__fasta_writer.write_multi(
						self.__fasta.get_seq(chrom, interval_chrom[0] + 1,
						interval_chrom[1]))

	def __write_chrom_full(self, chrom: int):
		"""Writes a chromosome to the opened Fasta file without breakpoints.
		:param chrom: Chromosome index
		"""
		self.__fasta_writer.set_bpl(
				self.__fasta.faidx.index[self.__fasta[chrom].name].lenc)
		self.__fasta_writer.write_header(self.__fasta[chrom].long_name)
		self.__fasta_writer.write_multi(self.__fasta[chrom])

	def __generate_all_breakpoints(self) -> dict[int, dict[str, list[int]]]:
		"""Generates all breakpoints for every chromosome pair.
		:return breakpoints: Dict containing all breakpoints for each
		chromosome and their partners
		"""
		breakpoints = {}

		for chrom in self.__get_pairs_to_mut_once():
			bp_amount = self.__get_bp_amount(
					len(self.__fasta[chrom]),
					self.__sim.chromosomes[chrom].it_rate,  #type:ignore
					len(self.__fasta[self.__partners[chrom]]),
					self.__sim.chromosomes[
					self.__partners[chrom]].it_rate)  #type:ignore

			bp_chrom, bp_partner = self.__get_breakpoints(
					chrom, len(self.__fasta[chrom]),
					len(self.__fasta[self.__partners[chrom]]), bp_amount)

			if (bp_chrom and bp_partner):
				breakpoints[chrom] = {"self": bp_chrom, "partner": bp_partner}
				breakpoints[self.__partners[chrom]] = {
						"self": bp_partner,
						"partner": bp_chrom
				}
			else:
				if not self.__args.ignore_warnings:
					print(f"{Colors.warn}WARNING: No interchromosomal translocations could be generated between sequence {chrom+1} and {self.__partners[chrom]+1} (it rates too low).{Colors.norm}",
							file=stderr)
		return breakpoints

	def __mutate_sequence(self, breakpoints: dict[int, dict[str, list[int]]]):
		"""Executes the interchromosomal translocations and writes them out.
		:param breakpoints: Dict of each chrom-partner breakpoint pairs
		"""
		for chrom in trange(len(self.__sim.chromosomes),
				desc="IT Mutating Sequences"):
			self.__fasta_writer.set_bpl(
					self.__fasta.faidx.index[self.__fasta[chrom].name].lenc)
			self.__fasta_writer.write_header(self.__fasta[chrom].long_name)

			if chrom in breakpoints:
				self.__write_with_bp(self.__fasta[chrom].name,
						breakpoints[chrom]["self"], len(self.__fasta[chrom]),
						self.__fasta[self.__partners[chrom]].name,
						breakpoints[chrom]["partner"],
						len(self.__fasta[self.__partners[chrom]]))
				self.__bedpe_writer.write(
						self.__fasta[chrom].name, breakpoints[chrom]["self"],
						len(self.__fasta[chrom]),
						self.__fasta[self.__partners[chrom]].name,
						breakpoints[chrom]["partner"],
						len(self.__fasta[self.__partners[chrom]]))
			else:
				self.__write_chrom_full(chrom)

	def mutate(self):
		"""Creates interchromosomal translocations using breakpoints
		and writes them to a Fasta and BEDPE file.
		"""
		breakpoints = self.__generate_all_breakpoints()
		self.__mutate_sequence(breakpoints)
