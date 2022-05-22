from __future__ import annotations

from hashlib import md5
from itertools import tee
from random import sample
from typing import TYPE_CHECKING, Any

from pyfaidx import Fasta

if TYPE_CHECKING:
	from pathlib import Path


class FastaDuplicateHeaderError(Exception):
	"""Raised when a Fasta contains at least one duplicate header."""


def get_md5(fname: Path) -> str:
	"""Calculates the md5 hash of a file.
	:param fname: Path to the file
	:return: MD5 hash
	"""
	hash_md5 = md5()
	with open(fname, "rb") as f:
		for chunk in iter(lambda: f.read(4096), b""):
			hash_md5.update(chunk)
	return hash_md5.hexdigest()


def load_fasta(fname: Path) -> Fasta:
	"""Loads a Fasta file as a pyfaidx.Fasta object.
	:param fname: Path to the file
	:return: Fasta object
	:raises FastaDuplicateHeaderError: When the Fasta contains duplicate header
	"""
	try:
		return Fasta(str(fname.resolve()),
				one_based_attributes=False,
				as_raw=True,
				sequence_always_upper=True,
				read_ahead=10000)
	except ValueError:
		raise FastaDuplicateHeaderError(
				f"Fasta {fname} contains duplicate header")


def sample_with_minimum_distance(n: int, k: int, d: int) -> list[int]:
	"""Sample of k elements from range(n), with a minimum distance d. Cannot be
	0 or n.
	:param n: Range 1-n
	:param k: Sample size
	:param d: Minimum distance
	:return: Random sample
	"""
	sampl = sample(range(1, n - (k - 1) * d), k)
	indices = sorted(range(len(sampl)), key=lambda i: sampl[i])
	return [
			s + d * r
			for s, r in zip(sampl, sorted(indices, key=lambda i: indices[i]))
	]


def pairwise(iterable: list[Any]) -> "zip"[Any]:
	"""Adjacent iteration.
	:param iterable: Iterable object
	:return: Adjacent values
	"""
	a, b = tee(iterable)
	next(b, None)
	return zip(a, b)
