from enum import Enum, auto


class MutType(Enum):
	"""Mutation types supported by Mutation-Simulator."""
	SN = auto()
	IN = auto()
	DE = auto()
	DU = auto()
	IV = auto()
	TL = auto()
	TLI = auto()
