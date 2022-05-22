from __future__ import annotations

# from functools import total_ordering
from random import randint, sample, shuffle, uniform
from sys import stderr
from typing import TYPE_CHECKING, Optional, Tuple

from numpy.random import choice  # type: ignore
from tqdm import tqdm

from .colors import Colors
from .fasta_writer import FastaWriter
from .mut_types import MutType
from .vcf_writer import VcfRecord, VcfWriter

if TYPE_CHECKING:
	from argparse import Namespace

	from pyfaidx import Fasta, FastaRecord

	from .rmt import RangeDefinition, SimulationSettings


# NOTE: Kept code if ordering is ever important
# @total_ordering
class Mutation:
	"""Holds information about a generated mutation."""
	__slots__ = ("type", "start", "stop", "trans_reverse", "trans_insert_pos")

	def __init__(self,
			type: MutType,
			start: int,
			stop: int = 0,
			trans_reverse: bool = False,
			trans_insert_pos: int = 0):
		"""Constructor. All positions are 0 based and inclusive.
		:param type: Mutation type
		:param start: Start position
		:param stop: Stop position
		:param trans_reverse: (Only necessary if type=TL) True if translocation is to be reversed
		:param trans_insert_pos: (Only necessary if type=TL) Position of the translocation insert
		"""
		self.type = type
		self.start = start
		self.stop = stop
		self.trans_reverse = trans_reverse
		self.trans_insert_pos = trans_insert_pos

	def __repr__(self) -> str:
		return f"{self.type}, {self.start}, {self.stop}, {self.trans_reverse}, {self.trans_insert_pos}\n"

	# def __le__(self, obj):
	# 	if self.type is MutType.TLI and obj.type is not MutType.TLI:
	# 		return ((self.trans_insert_pos) <= (obj.start))
	# 	elif self.type is MutType.TLI and obj.type is MutType.TLI:
	# 		return ((self.trans_insert_pos) <= (obj.trans_insert_pos))
	# 	elif self.type is not MutType.TLI and obj.type is MutType.TLI:
	# 		return ((self.start) <= (obj.trans_insert_pos))
	# 	else:
	# 		return ((self.start) <= (obj.start))

	# def __eq__(self, obj):
	# 	if self.type is MutType.TLI and obj.type is not MutType.TLI:
	# 		return (self.trans_insert_pos == obj.start)
	# 	elif self.type is MutType.TLI and obj.type is MutType.TLI:
	# 		return (self.trans_insert_pos == obj.trans_insert_pos)
	# 	elif self.type is not MutType.TLI and obj.type is MutType.TLI:
	# 		return (self.start == obj.trans_insert_pos)
	# 	else:
	# 		return (self.start == obj.start)


class Mutator:
	"""Stores required data for the simulation and performs mutations."""
	non_ambiguous = str.maketrans("KSYMWRBDHV-", "GCCAAACAAAN")
	complement = str.maketrans("ACGTUMRWSYKVHDB", "TGCAAKYWSRMBDHV")
	transitions = str.maketrans("AGTC", "GACT")

	def __init__(self, args: Namespace, fasta: Fasta, sim: SimulationSettings):
		"""Constructor.
		:param args: Commandline arguments
		:param fasta: Fasta file used in the simulation
		:param sim: Generated SimulationSettings from args, it or rmt mode
		"""
		self.__args = args
		self.__fasta = fasta
		self.__sim = sim

		self.__fasta_writer = FastaWriter(self.__args.outfasta)
		self.__vcf_writer = VcfWriter(self.__args.outvcf)
		self.__vcf_writer.write_header(self.__args.infile.name, self.__fasta,
				self.__sim.assembly_name, self.__sim.species_name,
				self.__sim.sample_name)

	def close(self):
		"""Closes the fasta_writer and vcf_writer filehandles. Does not close
		the fasta.
		"""
		self.__fasta_writer.close()
		self.__vcf_writer.close()

	def __del__(self):
		self.close()

	def mutate(self):
		"""Creates randon mutations and writes them to a Fasta and VCF file."""
		pbar = tqdm(total=len(self.__sim.chromosomes),
				desc="Mutating Sequences",
				position=0)
		for chrom in self.__sim.chromosomes:
			muts: dict[int, Mutation] = {}
			tls: list[int] = []
			tlis: list[int] = []

			for rng in chrom.range_definitions:
				if rng.mutation_settings.has_mutations:
					rng_muts, rng_tls, rng_tlis = self.__get_mutations(
							rng, len(self.__fasta[chrom.number]))
					if rng_muts:
						muts.update(rng_muts)
						tls.extend(rng_tls)
						tlis.extend(rng_tlis)

			if not muts and not self.__args.ignore_warnings:
				pbar.write(
						f"{Colors.warn}WARNING: No mutations could be generated on sequence {chrom.number+1} (mutation rates too low).{Colors.norm}",
						file=stderr)
			if tls:
				muts = self.__link_tls(muts, tls, tlis)

			self.__fasta_writer.set_bpl(self.__fasta.faidx.index[self.__fasta[
					chrom.number].name].lenc)
			self.__fasta_writer.write_header(
					self.__fasta[chrom.number].long_name)
			# Called even if no muts else the output fasta would miss
			# non mutated chromosomes
			self.__mutate_sequence(self.__fasta[chrom.number], muts,
					self.__sim.titv)
			pbar.update(1)
		pbar.close()

	def __get_mutations(
			self, rng: RangeDefinition, chrom_leng: int
	) -> Tuple[dict[int, Mutation], list[int], list[int]]:
		"""Generates a dictionary of position-mutation pairs. Translocations
		are additionally returned in separate lists.
		:param rng: RangeDefinition to simulate mutations on
		:param chrom_leng: Length of the sequence
		:return mutations: Dict with position-mutation pairs
		:return tl: List of translocational deletions
		:return tli: List of translocational insertions
		"""
		tls = []
		tlis = []
		start_positions = self.__get_mut_positions(
				rng.start,
				rng.stop,  #type: ignore
				sum(rng.mutation_settings.mut_rates.values()))  #type: ignore

		if not start_positions:
			return {}, [], []

		mutations = {
				start: Mutation(type=typ, start=start)
				for typ, start in tqdm(
				zip(
				choice(
				list(rng.mutation_settings.mut_chances.keys()),  # type: ignore
				p=list(
				rng.mutation_settings.mut_chances.values()),  # type: ignore
				size=len(start_positions)),
				start_positions),
				desc=
				f"Generating mutations in range: {rng.start+1}-{rng.stop+1}",  # type: ignore
				position=1,
				total=len(start_positions),
				leave=False)
		}

		last_mut_range = range(0)
		for pos in tqdm(list(mutations.keys()),
				desc="Assigning boundaries",
				position=1,
				leave=False):
			if pos in last_mut_range:
				del mutations[pos]
				continue

			mut = self.__get_stop_position(
					mutations[pos],
					rng.mutation_settings.mut_lengs,  # type: ignore
					chrom_leng)

			if not isinstance(mut, Mutation):
				del mutations[pos]
				continue

			mutations[pos] = mut
			if mut.type in [MutType.SN, MutType.IN]:
				last_mut_range = range(mut.start, mut.start + 1 +
						self.__sim.mut_block[mut.type])  #type: ignore
			else:
				last_mut_range = range(mut.start, mut.stop + 1 +
						self.__sim.mut_block[mut.type])  #type: ignore
				if mut.type is MutType.TL:
					tls.append(pos)
				if mut.type is MutType.TLI:
					tlis.append(pos)
		return mutations, tls, tlis

	@staticmethod
	def __get_mut_positions(start: int, stop: int,
			mut_rate: float) -> list[int]:
		"""Returns a sorted list of all mutation starting posititions.
		:param start: First base of the sequence
		:param stop: Last base of the sequence
		:param mut_rate: Mutation rate (sum)
		:return positions: List of mutation positions
		"""
		mut_count = int(((stop - start) + 1) * mut_rate)
		positions = sorted(sample(range(start, stop + 1), mut_count))
		return positions

	@staticmethod
	def __get_stop_position(mutation: Mutation, mut_lengs: dict[str,
			dict[MutType, int]], chrom_leng: int) -> Optional[Mutation]:
		"""Assigns a random stop value to a Mutation.
		:param mutation: Mutation without stop position
		:param mut_lengs: Mutation length settings
		:param chrom_leng: Length of the sequence
		:return mutation: Mutation with stop value or False if the mutation can
		not fit
		"""
		# Needs a minimum length of 2
		# If it cant fit at the end it will be skipped
		if mutation.type is MutType.SN:
			mutation.stop = mutation.start
		elif mutation.type is MutType.IV:
			if mutation.start + mut_lengs["max"][MutType.IV] >= chrom_leng - 1:
				return None
			else:
				mutation.stop = randint(
						mutation.start + mut_lengs["min"][MutType.IV] - 1,
						mutation.start + (mut_lengs["max"][MutType.IV] - 1))
		elif mutation.type is MutType.IN:
			mutation.stop = randint(
					mutation.start + mut_lengs["min"][mutation.type] - 1,
					mutation.start + mut_lengs["max"][mutation.type] - 1)
		elif mutation.type is MutType.DU:
			mutation.stop = randint(
					mutation.start + mut_lengs["min"][mutation.type] - 1,
					mutation.start + mut_lengs["max"][mutation.type] - 1)
			if mutation.stop > chrom_leng - 1:
				mutation.stop = chrom_leng - 1
		elif mutation.type in [MutType.TL, MutType.DE]:
			mutation.stop = randint(
					mutation.start + mut_lengs["min"][mutation.type] - 1,
					mutation.start + mut_lengs["max"][mutation.type] - 1)
			if mutation.stop > chrom_leng - 1:
				mutation.stop = chrom_leng - 1
		return mutation

	@classmethod
	def __link_tls(cls, muts: dict[int, Mutation], tls: list[int],
			tlis: list[int]) -> dict[int, Mutation]:
		"""Assigns every tl a random unique tli.
		:param muts: Dict of position-mutation pairs
		:param tls: List of translocational deletions
		:param tlis: List of translocational insertions
		:return muts: Dict of position-mutation pairs with tl/tli included and
		linked
		"""
		if len(tls) != len(tlis):
			tls, tlis = cls.__fix_tl_amount(tls, tlis)
		shuffle(tls)
		for tl_pos, tli_pos in zip(tls, tlis):
			muts[tli_pos] = Mutation(
					MutType.TLI, tl_pos, muts[tl_pos].stop,
					cls.__transloc_invert(muts[tl_pos].stop + 1 -
					muts[tl_pos].start), tli_pos)
		return muts

	@staticmethod
	def __fix_tl_amount(tls: list[int], tlis: list[int]):
		"""Removes tl or tli elements until both are equal in length.
		:param tls: List of translocational deletions
		:param tlis: List of translocational insertions
		:return tls: List of translocational deletions with equal length
		:return tlis: List of translocational insertions with equal length
		"""
		while len(tls) < len(tlis):
			del tlis[randint(0, len(tlis) - 1)]
		while len(tls) > len(tlis):
			del tls[randint(0, len(tls) - 1)]
		return tls, tlis

	@staticmethod
	def __transloc_invert(transloc_len: int) -> bool:
		"""Picks True or False randomly if the length of the sequence is >1.
		This is used to determine if a translocation will invert or not.
		:param transloc_len: Length of the translocation
		:return: True or False
		"""
		if randint(0, 1) == 0 or transloc_len < 2:
			return False
		else:
			return True

	def __mutate_sequence(self, sequence: FastaRecord, muts: dict[int,
			Mutation], titv: float):
		"""Mutates a sequence with a given mutation list and returns the new
		sequence and a list of VCF records.
		:param sequence: FastaRecord to mutate
		:param muts: Dict of position-mutation pairs
		:param titv: Transition/transversion rate
		"""
		pos = 0
		pbar = tqdm(total=len(sequence),
				desc="Writing",
				position=1,
				leave=False)
		while pos < len(sequence):
			if pos in muts:
				if muts[pos].type is MutType.SN:
					record = VcfRecord("sn",
							pos + 1,
							ref=self.__convert_ambiguous(sequence[pos]))
					record.alt = self.__get_snp(record.ref, titv)

					self.__fasta_writer.write(record.alt)
					self.__vcf_writer.write(record, sequence.name)

				elif muts[pos].type is MutType.IN:
					insert = self.__get_insert(muts[pos].stop + 1 - pos)
					record = VcfRecord("INS", pos, pos)
					if pos > 0:
						record.ref = self.__convert_ambiguous(sequence[pos -
								1])
						record.alt = f"{record.ref}{insert}"
					else:
						record.ref = self.__convert_ambiguous(sequence[0])
						record.alt = f"{insert}{record.ref}"
						record.start += 1
						record.end += 1
					record.len = len(insert)

					self.__fasta_writer.write_multi(f"{insert}{sequence[pos]}")
					self.__vcf_writer.write(record, sequence.name)

				elif muts[pos].type is MutType.DE or muts[
						pos].type is MutType.TL:
					svtype = "DEL" if muts[pos].type is MutType.DE else "DEL:ME"
					record = VcfRecord(svtype, pos, muts[pos].stop + 1)
					if pos > 0:
						record.ref = self.__convert_ambiguous(sequence[pos -
								1:record.end])
						record.alt = record.ref[0]
					else:
						record.start += 1
						record.end += 1
						record.ref = self.__convert_ambiguous(
								sequence[0:record.end])
						record.alt = record.ref[-1]
					record.len = muts[pos].stop - pos + 1

					pos = muts[pos].stop
					self.__vcf_writer.write(record, sequence.name)

				elif muts[pos].type is MutType.IV:
					record = VcfRecord("INV", pos + 1, muts[pos].stop + 1)
					record.ref = self.__convert_ambiguous(
							sequence[pos:record.end])
					record.alt = record.ref[::-1].translate(self.complement)

					self.__fasta_writer.write_multi(record.alt)
					pos = muts[pos].stop
					self.__vcf_writer.write(record, sequence.name)

				elif muts[pos].type is MutType.DU:
					record = VcfRecord("DUP", pos + 1)
					dupe = sequence[pos:muts[pos].stop + 1]
					record.ref = dupe
					record.alt = record.ref * 2
					record.len = len(dupe)
					record.end = pos + len(dupe)

					self.__fasta_writer.write_multi(record.alt)
					pos = muts[pos].stop
					self.__vcf_writer.write(record, sequence.name)

				elif muts[pos].type is MutType.TLI:
					record = VcfRecord("INS:ME")
					insert = self.__convert_ambiguous(
							sequence[muts[pos].start:muts[pos].stop + 1])
					if muts[pos].trans_reverse:
						insert = insert[::-1].translate(self.complement)
					record.start = pos
					if muts[pos].trans_insert_pos > 0:
						record.ref = self.__convert_ambiguous(sequence[pos -
								1:pos])
						record.alt = f"{record.ref}{insert}"
					else:
						record.start += 1
						record.ref = self.__convert_ambiguous(
								sequence[pos:pos + 1])
						record.alt = f"{insert}{record.ref}"
					record.end = record.start
					record.len = len(insert)

					self.__fasta_writer.write_multi(f"{insert}{sequence[pos]}")
					self.__vcf_writer.write(record, sequence.name)
			else:
				self.__fasta_writer.write(sequence[pos])
			pos += 1
			pbar.update(1)
		pbar.close()

	@classmethod
	def __get_snp(cls, base: str, ti_tv: float) -> str:
		"""Returns a mutated base with the probability given by the transition
		/ transversion rate.
		:param base: Nucleotide base
		:param titv: Transition/transversion rate
		:return: SNP
		"""
		p_ti = ti_tv * (1 / (ti_tv + 1))
		p = uniform(0, 1)
		if p <= p_ti:
			return cls.__get_ti_Base(base)
		else:
			return cls.__get_tv_Base(base)

	@staticmethod
	def __get_tv_Base(base: str) -> str:
		"""Returns one of the two transversion bases with equal probability.
		:param base: Nucleotide base
		:return: Transversion SNP
		"""
		return {
				"A": ["T", "C"],
				"G": ["C", "T"],
				"T": ["G", "A"],
				"C": ["A", "G"],
				"N": ["N", "N"]
		}[base][randint(0, 1)]

	@classmethod
	def __get_ti_Base(cls, base: str) -> str:
		"""Returns the transition base.
		:param base: Nucleotide base
		:return: Transition SNP
		"""
		return base.translate(cls.transitions)

	@staticmethod
	def __get_insert(leng: int) -> str:
		"""Returns a random sequence of specified length.
		:param leng: Length of the generated sequence
		:return: Random sequence of bases
		"""
		return "".join(choice(["A", "T", "G", "C"], leng))

	@classmethod
	def __convert_ambiguous(cls, substring: str) -> str:
		"""Returns the input string with converted ambiguity codes.
		:param substring: String to convert
		:return: String of bases without IUPAC codes
		"""
		return substring.translate(cls.non_ambiguous)
