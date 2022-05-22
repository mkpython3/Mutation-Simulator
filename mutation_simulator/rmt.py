from __future__ import annotations

from pathlib import Path
from sys import stderr
from typing import TYPE_CHECKING, Optional, Tuple

from .colors import Colors
from .defaults import Defaults
from .mut_types import MutType

if TYPE_CHECKING:
	from argparse import Namespace

	from pyfaidx import Fasta

MUT_INDICATORS = set(["sn", "in", "de", "iv", "du", "tl"])
MAX_LENG_INDICATORS = set(
		[f"{ind}max" for ind in MUT_INDICATORS if ind != "sn"])
MIN_LENG_INDICATORS = set(
		[f"{ind}min" for ind in MUT_INDICATORS if ind != "sn"])
META_KEYWORDS_STR = set([
		"fasta",
		"md5",
		"species_name",
		"assembly_name",
		"sample_name",
])
META_KEYWORDS_FLOAT = set(["titv"])
META_KEYWORDS_BLOCK = set([f"{ind}_block" for ind in MUT_INDICATORS])


class RatesTooHighError(Exception):
	"""Raised when the sum of all mutation rates is above 1."""


class RatesTooLowError(Exception):
	"""Raised when the sum of all mutation rates is 0 or below."""


class ItRateTooHighError(Exception):
	"""Raised when the specified it rate is higher than 1."""


class ItRateTooLowError(Exception):
	"""Raised when the specified it rate is below 0."""


class ITNotEnoughAvailChromsError(Exception):
	"""Raised when not enought chromosomes are available."""


class TitvTooLowError(Exception):
	"""Raised when the specified titv ratio is below 0."""


class RMTParseError(Exception):
	"""Raised when the sum of all mutation rates is above 1."""


class MissingLengthError(Exception):
	"""Raised when a mutation rate is specified without a min or max length."""


class MinimumLengthTooLowError(Exception):
	"""Raised when the specified minimum length for a mut_type is too low."""


class MinimumLengthHigherThanMaximumError(Exception):
	"""Raised when the min length is higher than the max length or vice versa."""


class ChromNotExistError(Exception):
	"""Raised when a chromosome index in the RMT file does not exist in the fasta file."""


class RangeDefinitionOutOfBoundsError(Exception):
	"""Raised when a range definition is out of bounds of the chromosome start/end values."""


class MutationSettings:
	"""Holds details of mutation rate and length settings."""

	def __init__(self, mut_rates: Optional[dict[MutType, float]],
			mut_lengs: Optional[dict[str, dict[MutType, int]]]):
		"""Constructor.
		:param mut_rates: Mutation rates per type or None
		:param mut_lengs: Mutation min and max lengths per type or None
		"""
		self.mut_rates = mut_rates
		self.mut_lengs = mut_lengs
		self.__validate()
		if self.mut_rates:
			if MutType.TL in self.mut_rates:
				self.mut_rates[MutType.TL] = self.mut_rates[MutType.TL] / 2
				self.mut_rates[MutType.TLI] = self.mut_rates[MutType.TL]
		self.__mut_rates_2_chances()

	def __validate(self):
		"""Checks if specified rate and length settings are valid.
		:raises RatesTooHighError: If a / all rate/s are too high
		:raises RatesTooLowError: If a / all rate/s are too low
		:raises MissingLengthError: If a rate is specified without min or max length (except sn)
		:raises MinimumLengthHigherThanMaximumError: If a minimum length is higher than its corresponding maximum length
		:raises MinimumLengthTooLowError: If a minimum length is lower than possible
		"""
		if self.mut_rates:
			if any(rate < 0 for rate in self.mut_rates.values()) or sum(
					self.mut_rates.values()) <= 0:
				raise RatesTooLowError(
						f"Mutation rate/s too low for: {', '.join(typ.name for typ,rate in self.mut_rates.items() if rate <=0)}. If this was intentional use the None keyword instead"
				)
			if sum(self.mut_rates.values()) > 1:
				too_high_rates = ", ".join(typ.name
						for typ, rate in self.mut_rates.items() if rate > 1)
				if too_high_rates:
					raise RatesTooHighError(
							f"Mutation rate/s too high for: {too_high_rates}")
				else:
					raise RatesTooHighError("Sum of mutation rates too high")
			for mut_type in self.mut_rates:
				if mut_type is not MutType.SN:
					if not self.mut_lengs or mut_type not in self.mut_lengs[
							"min"] or mut_type not in self.mut_lengs["max"]:
						raise MissingLengthError(
								f"Missing length keyword for: {mut_type.name}")
					if self.mut_lengs["min"][mut_type] > self.mut_lengs["max"][
							mut_type]:
						raise MinimumLengthHigherThanMaximumError(
								f"Minimum length is greater and maximum length for: {mut_type.name}"
						)
					if mut_type is not MutType.IV and self.mut_lengs["min"][
							mut_type] < 1:
						raise MinimumLengthTooLowError(
								f"Minimum length too low for: {mut_type.name}")
					if mut_type is MutType.IV and self.mut_lengs["min"][
							mut_type] < 2:
						raise MinimumLengthTooLowError(
								f"Minimum length too low for: {mut_type.name}")

	def __repr__(self) -> str:
		return f"Rates: {self.mut_rates}, Chances: {self.mut_chances}, Lengs: {self.mut_lengs}"

	def __mut_rates_2_chances(self):
		"""Calculates the chances of each mutation as if there would be only
		one mutation rate.
		"""
		if self.mut_rates:
			rate_sum = sum(self.mut_rates.values())
			self.mut_chances = {
					mut_type: chance / rate_sum
					for mut_type, chance in self.mut_rates.items()
			}
		else:
			self.mut_chances = None

	@property
	def has_mutations(self) -> bool:
		"""Checks if any mutation rate is higher than 0.
		:return: True if any mutation rate is non zero else False
		"""
		if self.mut_rates and any(rates for rates in self.mut_rates.values()):
			return True
		else:
			return False


class RangeDefinition:
	"""Specifies the range in which the MutationSettings should apply."""

	def __init__(self, start: int, stop: int | str,
			mutation_settings: MutationSettings):
		"""Constructor.
		:param start: Start of the range (0 based, inclusive)
		:param stop: Stop of the range (0 based, inclusive)
		:param mutation_settings: MutationSettings object
		"""
		self.start = start
		self.stop = stop
		self.mutation_settings = mutation_settings

	def __repr__(self) -> str:
		return f"{self.start}-{self.stop} {self.mutation_settings}"

	def eval_end(self, chrom_length: int):
		"""Replaces the 'end' stop attribute value with the specified chromosome length.
		:param chrom_length: Length of a chromosome
		:notes: Do not input length-1
		"""
		if self.stop == "end":
			self.stop = chrom_length - 1


class ChromosomeSettings:
	"""Holds the it rate and range definitions of each chromosome."""

	def __init__(self, number: int, it_rate: Optional[float],
			range_definitions: list[RangeDefinition]):
		"""Constructor.
		:param number: Index of the chromosome in the Fasta
		:param it_rate: Interchromosomal translocation rate of None
		:param range_definitions: List of RangeDefinition objects
		"""
		self.number = number
		self.it_rate = it_rate
		self.__validate_it()
		self.range_definitions = range_definitions

	def __repr__(self) -> str:
		return f"{self.number}\nit={self.it_rate}\n{self.range_definitions}\n"

	def __validate_it(self):
		"""Checks if the specified it rate is valid.
		:raises ItRateTooLowError: If the rate is too low
		:raises ItRateTooHighError: If the rate is too high
		"""
		if self.it_rate and self.it_rate > 0.5:
			raise ItRateTooHighError(
					"Interchromosomal translocation rate too high")
		if self.it_rate and self.it_rate < 0:
			raise ItRateTooLowError(
					"Interchromosomal translocation rate too low")

	def fill_missing_ranges(
			self, chrom_length: int, std_mut_settings: MutationSettings):
		"""Fills range definitions that are not specified with the standard.
		:param chrom_length: Length of the chromosome
		:param std_mut_settings: Standard MutationSettings
		"""
		if self.range_definitions:
			# Check if starts at 0
			if self.range_definitions[0].start != 0:
				self.range_definitions.insert(
						0,
						RangeDefinition(0, self.range_definitions[0].start - 1,
						std_mut_settings))
			# Check if it ends at chrom_length-1
			if self.range_definitions[-1].stop != chrom_length - 1:
				self.range_definitions.append(
						RangeDefinition(
						self.range_definitions[-1].stop + 1,  # type: ignore
						chrom_length - 1,
						std_mut_settings))
			# Check inbetween ranges
			if len(self.range_definitions) > 1:
				i = 0
				while i < len(self.range_definitions) - 1:
					if not self.range_definitions[  # type: ignore
							i].stop + 1 == self.range_definitions[i + 1].start:
						self.range_definitions.insert(
								i + 1,
								RangeDefinition(
								self.range_definitions[i].stop +
								1,  # type: ignore
								self.range_definitions[i + 1].start - 1,
								std_mut_settings))
					i += 1
		else:
			self.range_definitions.append(
					RangeDefinition(0, chrom_length - 1, std_mut_settings))


class SimulationSettings:
	"""Holds every setting for a mutation simulation."""

	def __init__(self,
			std: MutationSettings,
			std_it: Optional[float],
			chromosomes: list[ChromosomeSettings],
			mut_block: Optional[dict[MutType, int]],
			fasta: Optional[Path] = None,
			md5: Optional[str] = None,
			titv: float = Defaults.TITV,
			species_name: str = Defaults.SPECIES_NAME,
			assembly_name: str = Defaults.ASSEMBLY_NAME,
			sample_name: str = Defaults.SAMPLE_NAME,
			ignore_warnings: bool = Defaults.IGNORE_WARNINGS):
		"""Constructor.
		:param std: Standard MutationSettings
		:param std_it: Standard interchromosomal translocation rate or None
		:param chromosomes: List of ChromosomeSettings
		:param mut_block: Block values per type or None
		:param fasta: Path to the Fasta file or None
		:param md5: MD5 hash of the Fasta file or None
		:param titv: Transition/transversion rate
		:param species_name: Name of the species
		:param assembly_name: Name of the assembly
		:param sample_name: Name of the sample
		:param ignore_warnings: Supresses warnings
		"""
		self.__std = std
		self.__std_it = std_it
		self.chromosomes = chromosomes
		self.mut_block = mut_block
		self.__validate_mut_block(ignore_warnings)
		self.fasta: Optional[Path] = fasta
		self.md5: Optional[str] = md5
		self.titv: float = titv
		self.__validate_titv()
		self.species_name: str = species_name
		self.assembly_name: str = assembly_name
		self.sample_name: str = sample_name

	def __repr__(self) -> str:
		return f"[META]\nfasta={self.fasta}\nmd5={self.md5}\ntitv={self.titv}\nspecies_name={self.species_name}\nassembly_name={self.assembly_name}\nsample_name={self.sample_name}\nmut_block={self.mut_block}\n\n[STD]\nit={self.__std_it}\n{self.__std}\n\n[RD]\n{self.chromosomes}"

	def __validate_it(self, fasta: Fasta):
		"""Checks if the sum of all it rates is above 0 and for enough
		available chromosomes.
		:param fasta: Fasta file used in the simulation
		"""
		avail_chroms = [
				chrom.number for chrom in self.chromosomes if
				chrom.it_rate is not None  # intentional is comparison, 0 is ok
				and len(fasta[chrom.number]) > 2
		]
		if len(avail_chroms) < 2:
			raise ITNotEnoughAvailChromsError(
					"Not enought available chromosomes for interchromosomal translocations"
			)
		if sum(self.chromosomes[n].it_rate
				for n in avail_chroms) == 0:  #type:ignore
			raise ItRateTooLowError(
					"Interchromosomal translocation rates are too low")

	def __validate_mut_block(self, ignore_warnings: bool):
		"""Checks if the mutation block values are valid. Adjusts to 1 if lower.
		:param ignore_warnings: Supresses warnings
		"""
		if self.mut_block:
			for mut_type in MutType:
				if mut_type is not MutType.TLI:
					if mut_type not in self.mut_block:
						self.mut_block[mut_type] = 1
					elif self.mut_block[mut_type] < 1:
						self.mut_block[mut_type] = 1
						if not ignore_warnings:
							print(f"{Colors.warn}WARNING: '{mut_type.name}' block was set to 1{Colors.norm}",
									file=stderr)
			# tli will never be set by user
			self.mut_block[MutType.TLI] = self.mut_block[MutType.TL]
		else:
			self.mut_block = Defaults.MUT_BLOCK

	def __validate_titv(self):
		"""Checks if the titv rate is valid.
		:raises TitvTooLowError: If the titv rate is too low
		"""
		if self.titv < 0:
			raise TitvTooLowError("Titv value is below 0")

	def __check_chroms_exist(self, fasta: Fasta):
		"""Checks if all SimulationSettings chromosomes exist in the Fasta.
		:param fasta: Fasta file used in the simulation
		:raises ChromNotExistError: If a chromosome does not exist
		"""
		existing_chroms = list(fasta.keys())
		for chrom in self.chromosomes:
			try:
				existing_chroms[chrom.number]
			except IndexError:
				raise ChromNotExistError(
						f"Chromosome {chrom.number+1} does not exist in the fasta file"
				)

	def __eval_chrom_ends(self, fasta: Fasta):
		"""Replaces the 'end' stop attribute values with the corresponding chromosome lengths.
		:param fasta: Fasta file used in the simulation
		"""
		for chrom in self.chromosomes:
			if chrom.range_definitions:
				chrom.range_definitions[-1].eval_end(len(fasta[chrom.number]))

	def __fill_missing_chroms(self, fasta: Fasta):
		"""Fills chromosomes not contained yet with standard values.
		:param fasta: Fasta file used in the simulation
		"""
		all_chroms = list(fasta.keys())
		self_chroms = set([chrom.number for chrom in self.chromosomes])
		missing_chroms = [
				chrom_idx for chrom_idx, _ in enumerate(all_chroms)
				if chrom_idx not in self_chroms
		]
		if missing_chroms:
			for chrom_idx in missing_chroms:
				self.chromosomes.append(
						ChromosomeSettings(chrom_idx, self.__std_it, [
						RangeDefinition(0,
						len(fasta[chrom_idx]) - 1, self.__std)
						]))

	def __sort_range_definitions(self):
		"""Sorts chromosomes and their range definitions."""
		self.chromosomes = sorted(self.chromosomes,
				key=lambda chrom: chrom.number)
		for i, _ in enumerate(self.chromosomes):
			self.chromosomes[i].range_definitions = sorted(
					self.chromosomes[i].range_definitions,
					key=lambda rd: rd.start)

	def __check_range_definitions_in_bounds(self, fasta):
		"""Checks if the boundaries of the first and last range definition are valid.
		:param fasta: Fasta file used in the simulation
		:raises RangeDefinitionOutOfBoundsError: If a range definition is < 0 or > chromosome length
		"""
		existing_chroms = list(fasta.keys())
		for chrom in self.chromosomes:
			if chrom.range_definitions:
				if chrom.range_definitions[0].start < 0:
					raise RangeDefinitionOutOfBoundsError(
							f"A range definition of chromosome {chrom.number+1} is starting at 0"
					)
				if chrom.range_definitions[-1].stop > len(
						fasta[existing_chroms[chrom.number]]):  #type: ignore
					raise RangeDefinitionOutOfBoundsError(
							f"A range definition of chromosome {chrom.number+1} is longer than the chromosome"
					)

	def __fill_missing_chrom_ranges(self, fasta: Fasta):
		"""Fills missing range definitions of all chromosomes.
		:param fasta: Fasta file used in the simulation
		"""
		for chrom in self.chromosomes:
			chrom.fill_missing_ranges(len(fasta[chrom.number]), self.__std)

	@classmethod
	def from_args(cls, args: Namespace, fasta: Fasta,
			ignore_warnings: bool) -> "SimulationSettings":
		"""Alternative constructor for the args mode.
		:param args: Commandline arguments
		:param fasta: Fasta file used in the simulation
		:param ignore_warnings: Supresses warnings
		:return: SimulationSettings with all chromosomes according to the commandline arguments
		:raises RatesTooHighError: If a / all rate/s are too high
		:raises RatesTooLowError: If a / all rate/s are too low
		:raises MissingLengthError: If a rate is specified without min or max length (except sn)
		:raises MinimumLengthHigherThanMaximumError: If a minimum length is higher than its corresponding maximum length
		:raises MinimumLengthTooLowError: If a minimum length is lower than possible
		:raises TitvTooLowError: If the titv rate is too low
		"""
		mut_rates = {
				MutType.SN: args.snp,
				MutType.IN: args.insert,
				MutType.DE: args.deletion,
				MutType.IV: args.inversion,
				MutType.DU: args.duplication,
				MutType.TL: args.translocation
		}
		mut_lengs = {
				"min": {
				MutType.IN: args.insertminlength,
				MutType.DE: args.deletionminlength,
				MutType.IV: args.inversionminlength,
				MutType.DU: args.duplicationminlength,
				MutType.TL: args.translocationminlength
				},
				"max": {
				MutType.IN: args.insertmaxlength,
				MutType.DE: args.deletionmaxlength,
				MutType.IV: args.inversionmaxlength,
				MutType.DU: args.duplicationmaxlength,
				MutType.TL: args.translocationmaxlength
				}
		}
		mut_block = {
				MutType.SN: args.snpblock,
				MutType.IN: args.insertblock,
				MutType.DE: args.deletionblock,
				MutType.IV: args.inversionblock,
				MutType.DU: args.duplicationblock,
				MutType.TL: args.translocationblock
		}
		sim = cls(MutationSettings(mut_rates, mut_lengs),
				None, [],
				mut_block,
				titv=args.transitionstransversions,
				species_name=args.species,
				assembly_name=args.assembly,
				sample_name=args.sample,
				ignore_warnings=ignore_warnings)
		sim.__fill_missing_chroms(fasta)
		sim.__sort_range_definitions()
		return sim

	@classmethod
	def from_it(cls, it_rate: float, fasta: Fasta,
			ignore_warnings: bool) -> "SimulationSettings":
		"""Alternative constructor for the it mode.
		:param it_rate: Interchromosomal translocation rate
		:param fasta: Fasta file used in the simulation
		:param ignore_warnings: Supresses warnings
		:return: SimulationSettings with all chromosomes having an interchromosomal translocation rate of it_rate
		:raises ItRateTooLowError: If the rate is too low
		:raises ItRateTooHighError: If the rate is too high
		:raises RangeDefinitionOutOfBoundsError: If a range definition is < 0 or > chromosome length
		"""
		sim = cls(MutationSettings(None, None),
				it_rate, [],
				None,
				ignore_warnings=ignore_warnings)
		sim.__fill_missing_chroms(fasta)
		sim.__sort_range_definitions()
		sim.__validate_it(fasta)
		return sim

	@classmethod
	def from_rmt(cls, path: str | Path, fasta: Fasta,
			ignore_warnings: bool) -> "SimulationSettings":
		"""Alternative constructor for the rmt mode.
		:param path: Path to the RMT file
		:param fasta: Fasta file used in the simulation
		:param ignore_warnings: Supresses warnings
		:return: SimulationSettings with all chromosomes and settings according to the RMT file
		:raises RatesTooHighError: If a / all rate/s are too high
		:raises RatesTooLowError: If a / all rate/s are too low
		:raises MissingLengthError: If a rate is specified without min or max length (except sn)
		:raises MinimumLengthHigherThanMaximumError: If a minimum length is higher than its corresponding maximum length
		:raises MinimumLengthTooLowError: If a minimum length is lower than possible
		:raises ItRateTooLowError: If the rate is too low
		:raises ItRateTooHighError: If the rate is too high
		:raises TitvTooLowError: If the titv rate is too low
		:raises ChromNotExistError: If a chromosome does not exist
		:raises RangeDefinitionOutOfBoundsError: If a range definition is < 0 or > chromosome length
		:raises RMTParseError: When a keyword value was expected to be numeric but is not interpretable
		"""
		rmt_raw = cls.__read_rmt_raw(path)
		rmt_preprocessed = cls.__split_raw_rmt(rmt_raw)
		if len(rmt_preprocessed["std"]) != 2:
			raise RMTParseError(
					f"Standard section not defined or malformed. Occurred while reading {path}"
			)
		try:
			meta, mut_block = cls.__parse_meta(rmt_preprocessed["meta"])
			std_it, std_mut_settings = cls.__parse_std(rmt_preprocessed["std"])
			chromosome_settings = cls.__parse_rd(rmt_preprocessed["rd"],
					std_it)

			sim = cls(std_mut_settings,
					std_it,
					chromosome_settings,
					mut_block,
					ignore_warnings=ignore_warnings,
					**meta)
			sim.__check_chroms_exist(fasta)
			sim.__eval_chrom_ends(fasta)
			sim.__fill_missing_chroms(fasta)
			sim.__sort_range_definitions()
			sim.__check_range_definitions_in_bounds(fasta)
			sim.__fill_missing_chrom_ranges(fasta)
			if sim.has_it:
				sim.__validate_it(fasta)
		except Exception as e:
			raise type(e)(f"{e}. Occurred while reading {path}")

		return sim

	@staticmethod
	def __read_rmt_raw(rmt_path: str | Path) -> list[str]:
		"""Reads an RMT file.
		:param rmt_path: Path to the RMT file
		:return rmt_raw: List of RMT lines without comments or empty lines
		"""
		rmt_raw = []
		with open(rmt_path, "r") as hndl:
			for line in hndl.readlines():
				if not line.startswith("#"):
					line = line.split("#")[0].strip()
					if line:
						rmt_raw.append(line)
		return rmt_raw

	@staticmethod
	def __split_raw_rmt(rmt_raw: list[str]) -> dict[str, list[str]]:
		"""Splits the rmt_raw into the RMT sections.
		:param rmt_raw: List of RMT lines without comments or empty lines
		:return rmt_preprocessed: Dict of RMT sections
		"""
		rmt_preprocessed = {"meta": [], "std": [], "rd": []}
		curr_section = "meta"

		for line in rmt_raw:
			line = line.lower()
			if line == "std":
				curr_section = "std"
				continue
			elif line.startswith("chr"):
				curr_section = "rd"

			rmt_preprocessed[curr_section].append(line)
		return rmt_preprocessed

	@staticmethod
	def __parse_meta(
			rmt_meta: list[str]
	) -> Tuple[dict[str, str | float], dict[MutType, int]]:
		"""Evaluates Meta section keywords.
		:param rmt_meta: RMT lines of the Meta section
		:return meta: Keyword value dict of the specified Metakeys
		:return mut_block: Separate mutation block values specified
		:raises RMTParseError: When a keyword value was expected to be numeric but is not interpretable
		"""
		meta = {}
		mut_block = {}

		for line in rmt_meta:
			key, val = [tok.strip() for tok in line.split("=")]
			if key in META_KEYWORDS_STR:
				meta[key] = val
			elif key in META_KEYWORDS_FLOAT:
				try:
					meta[key] = float(val)
				except ValueError:
					raise RMTParseError(
							f"{key.capitalize()} value of '{val}' is not representable as a float"
					)
			elif key in META_KEYWORDS_BLOCK:
				mut_type = MutType[key.removesuffix("_block").upper()]
				try:
					mut_block[mut_type] = int(val)
				except ValueError:
					raise RMTParseError(
							f"Mut block value of {key} is not representable as an integer"
					)
		return meta, mut_block

	@classmethod
	def __parse_std(
			cls,
			rmt_std: list[str]) -> Tuple[Optional[float], MutationSettings]:
		"""Parses interchromosomal translocation rate and mutation settings.
		:param rmt_std: RMT lines of the Standard section
		:return it: Interchromosomal translocation rate of the Standard
		:return MutationSettings: Mutation rate and length settings of the Standard
		"""
		return cls.__parse_it(rmt_std[0]), cls.__parse_settings(rmt_std[1])

	@staticmethod
	def __parse_it(line: str) -> Optional[float]:
		"""Parses the interchromosomal translocation rate in a RMT line.
		:param line: RMT line
		:return it: Interchromosomal translocation rate
		:raises RMTParseError: When the format of the line is not according to the RMT spec.
		"""
		it = None

		tokens = [tok.strip() for tok in line.split(" ") if tok.strip()]
		if tokens[0] != "it" or len(tokens) != 2:
			raise RMTParseError(
					"Malformed interchromosomal translocation rate setting")
		else:
			if tokens[1] != "none":
				try:
					it = float(tokens[1])
				except ValueError:
					raise RMTParseError(
							"Malformed interchromosomal translocation rate setting"
					)
		return it

	@staticmethod
	def __parse_settings(line: str) -> MutationSettings:
		"""Parses the MutationSettings of a RMT line.
		:param line: RMT line
		:return MutationSettings: Mutation rate and length settings
		:raises RMTParseError: When the format of the line is not according to the RMT spec.
		"""
		if line == "none":
			return MutationSettings(None, None)

		tokens = [tok.strip() for tok in line.split(" ") if tok.strip()]
		if len(tokens) % 2 != 0:
			raise RMTParseError("Malformed mutation settings")

		mut_rates = {}
		mut_lengs = {"min": {}, "max": {}}

		for i, tok in enumerate(tokens[:-1]):
			try:
				if tok in MUT_INDICATORS:
					mut_rates[MutType[tok.upper()]] = float(tokens[i + 1])
				elif tok in MAX_LENG_INDICATORS:
					mut_lengs["max"][MutType[tok.removesuffix(
							"max").upper()]] = int(tokens[i + 1])
				elif tok in MIN_LENG_INDICATORS:
					mut_lengs["min"][MutType[tok.removesuffix(
							"min").upper()]] = int(tokens[i + 1])
			except ValueError:
				raise RMTParseError("Malformed mutation settings")

		return MutationSettings(mut_rates, mut_lengs)

	@staticmethod
	def __parse_range(range_str: str) -> Tuple[int, int | str]:
		"""Parses the range of a RMT line.
		:param range_str: Start-Stop string of a RMT line
		:return start: Start of the range
		:return stop: Stop of the range
		:raises RMTParseError: When the format of the line is not according to the RMT spec.
		"""
		tokens = [tok.strip() for tok in range_str.split("-") if tok.strip()]
		if len(tokens) != 2 or range_str.count("-") != 1:
			raise RMTParseError("Malformed range in range definitions")
		try:
			start = int(tokens[0]) - 1
			stop = tokens[1]
			if tokens[1] != "end":
				stop = int(tokens[1]) - 1
		except ValueError:
			raise RMTParseError("Malformed range in range definitions")
		return start, stop

	@classmethod
	def __parse_rd(cls, rmt_rd: list[str],
			std_it: Optional[float]) -> list[ChromosomeSettings]:
		"""Parses the Range definitions section.
		:param rmt_std: RMT lines of the Range definitions section
		:param std_it: Interchromosomal translocation rate of the Standard section
		:return start: Chromosome settings of the Range definitions section
		:raises RMTParseError: When the format of the line is not according to the RMT spec.
		"""
		# This dict is only used for parsing
		chromosomes = {}
		chrom_index = None

		for row in rmt_rd:
			if row.startswith("chr"):
				try:
					chrom_index = int(row.split(" ")[-1]) - 1
				except ValueError:
					raise RMTParseError(f"Chromosome index {row} is invalid")
				chromosomes[chrom_index] = {"range_definitions": []}
			elif row.startswith("it"):
				chromosomes[chrom_index]["it"] = cls.__parse_it(row)
			else:
				_range = row.partition(" ")[0]
				start, stop = cls.__parse_range(_range)
				settings = row.removeprefix(_range + " ")
				mut_settings = cls.__parse_settings(settings)
				chromosomes[chrom_index]["range_definitions"].append(
						RangeDefinition(start, stop, mut_settings))

		chromosome_settings = []
		for chrom in chromosomes:
			# Add missing std_it
			if "it" not in chromosomes[chrom]:
				chromosomes[chrom]["it"] = std_it

			# Construct ChromosomeSettings class from dict entry and add to return list
			chromosome_settings.append(
					ChromosomeSettings(chrom, chromosomes[chrom]["it"],
					chromosomes[chrom]["range_definitions"]))

		return chromosome_settings

	@property
	def has_mutations(self) -> bool:
		"""Checks if any mutation rate is higher than 0.
		:return: True if any mutation rate is non zero else False
		"""
		for chrom in self.chromosomes:
			if any(rd.mutation_settings.has_mutations
					for rd in chrom.range_definitions):
				return True
		return False

	@property
	def has_it(self) -> bool:
		"""Checks if any interchromosomal translocation rate is higher than 0.
		:return: True if any interchromosomal translocation rate is non zero else False
		"""
		if any(chrom.it_rate for chrom in self.chromosomes):
			return True
		else:
			return False
