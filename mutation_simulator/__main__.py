#!/usr/bin/env python3

# Mutation-Simulator
# Copyright (C) 2024 Marius Kühl

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

from __future__ import annotations

from timeit import default_timer as timer
from typing import TYPE_CHECKING, Tuple

from pyfaidx import FastaIndexingError, FastaNotFoundError

from . import *

if TYPE_CHECKING:
	from argparse import Namespace

	from pyfaidx import Fasta


def initialize() -> Tuple[Namespace, Fasta, SimulationSettings]:
	"""Initializes needed objects for the simulation and exits on error.
	:return args: Commandline argument object
	:return fasta: Genomic sequence
	:return sim: Settings object for the simulation
	"""
	args = get_args()
	try:
		fasta = load_fasta(args.infile)
		if args.mode == "args":
			sim = SimulationSettings.from_args(args, fasta,
					args.ignore_warnings)
		elif args.mode == "it":
			sim = SimulationSettings.from_it(args.interchromosomalrate, fasta,
					args.ignore_warnings)
		else:
			sim = SimulationSettings.from_rmt(args.rmtfile, fasta,
					args.ignore_warnings)
	except (FileNotFoundError, ITNotEnoughAvailChromsError, RatesTooHighError,
			RatesTooLowError, FastaIndexingError, FastaNotFoundError,
			ItRateTooHighError, ItRateTooLowError, RMTParseError,
			MissingLengthError, MinimumLengthTooLowError, TitvTooLowError,
			ChromNotExistError, RangeDefinitionOutOfBoundsError,
			FastaDuplicateHeaderError,
			MinimumLengthHigherThanMaximumError) as e:
		exit_with_error(e, args.no_color)
	if not args.ignore_warnings:
		warn_user(args, sim)
	return args, fasta, sim


def warn_user(args: Namespace, sim: SimulationSettings):
	"""Warns the user if the potentially wrong RMT or Fasta is used.
	:param args: Commandline arguments
	:param sim: Generated SimulationSettings from args, it or rmt mode
	"""
	if sim.fasta and args.infile.name != sim.fasta:
		print_warning("Fasta filename does not match RMT", args.no_color)
	if sim.md5 and get_md5(args.infile) != sim.md5:
		print_warning("Fasta md5 hash does not match RMT", args.no_color)


def main():
	"""Called when Mutation-Simulator is run."""
	start = timer()
	args, fasta, sim = initialize()
	if sim.has_mutations:
		try:
			mutator = Mutator(args, fasta, sim)
			mutator.mutate()
			mutator.close()
			fasta.close()
		except (FastaWriterError, VcfWriterError) as e:
			exit_with_error(e, args.no_color)
	if sim.has_it:
		# Reload after mutations
		if sim.has_mutations:
			try:
				fasta = load_fasta(args.outfasta)
			except (FastaDuplicateHeaderError, FastaIndexingError,
					FastaNotFoundError) as e:
				exit_with_error(e, args.no_color)
		try:
			it_mutator = ITMutator(args, fasta, sim)
			it_mutator.mutate()
			it_mutator.close()
			fasta.close()
		except (FastaWriterError, BedpeWriterError) as e:
			exit_with_error(e, args.no_color)
	stop = timer()
	runtime = round(stop - start, 4)
	if not args.quiet:
		print_success(f"Mutation-Simulator finished in: {runtime}s",
				args.no_color)


if __name__ == "__main__":
	main()
