from sys import stderr, stdout


class Colors:
	"""Holds color codes when in TTY mode."""
	if stdout.isatty() and stderr.isatty():
		error = "\033[0;31m"
		warn = "\033[0;33m"
		ok = "\033[0;32m"
		norm = "\033[0m"
	else:
		error = ""
		warn = ""
		ok = ""
		norm = ""
