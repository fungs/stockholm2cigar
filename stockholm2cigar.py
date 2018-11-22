#!/usr/bin/env python3
"""
Convert Stockholm format alignments as produced by HMMER v3.1b hmmalign to simple
CIGAR strings for each input sequences so that the columns of the original profile
hidden markov model alignment can be reconstructed.

Note: insertions are represented using soft masking
      (lower case symbols and * in alignment notation)
"""

__author__ = "code@fungs.de"
__version__ = "0.1.0"
__license__ = "GPLv3"

from collections import OrderedDict  # keeps order of input seqs
from copy import deepcopy

def read_stockholm(infile):
	"Parse Stockholm multiple sequence alignment format"
	
	try:
		line = infile.readline().rstrip("\n")
		assert (line == "# STOCKHOLM 1.0"), "File does not comply with the Stockholm format specs."
	except StopIteration:  # empty files not allowed
		raise IOError("Input is empty.")
	
	data = OrderedDict()
	
	for line in infile:
		line = line.rstrip()
		
		if not line:  # skip empty lines
			continue
		
		if line == '//':  # end of alignment
			return data
		
		if not line.startswith('#'):  # skip comments
			ident, seq = line.split(maxsplit=2)[:2]
			try:
				data[ident].extend(seq)
			except KeyError:
				data[ident] = list(seq)
	
	raise IOError("End of alignment could not be found.")


class CigarString(object):
	"""Class to construct and represent CIGAR strings"""
	
	def __init__(self):
		self._duplication = []
		self._symbol = []
		self._len = 0
	
	def append(self, s):
		if not self._len or s != self._symbol[-1]:
			self._symbol.append(s)
			self._duplication.append(1)
		else:
			self._duplication[-1] += 1
		
		self._len += 1
		return self
	
	def clear(self):
		self._len = 0
		del self._duplication[:]
		del self._symbol[:]
	
	def join(self, cigarstr):
		if len(cigarstr) == 0:
			return self
		if not self._len:
			self = deepcopy(cigarstr)
			return self
		
		self._len += len(cigarstr)
		
		if self._symbol[-1] == cigarstr._symbol[0]:
			self._duplication[-1] += cigarstr._duplication[0]
			self._symbol.extend(cigarstr._symbol[1:])
			self._duplication.extend(cigarstr._duplication[1:])
			return self
		
		self._symbol.extend(cigarstr._symbol)
		self._duplication.extend(cigarstr._duplication)
		return self
	
	def __add__(self, cigarstr):
		return deepcopy(self).join(cigarstr)
	
	@property
	def code(self):
		return "".join("{}{}".format(i, s) for i, s in zip(self._duplication, self._symbol))
	
	__len__ = lambda self: self._len
	__repr__ = lambda self: self.code


def cigar_from_aligned_hmmalign(alnseq):
	"""Take a sequence with SAM notation symbols (-.) and build a cigar string."""
	
	cs = CigarString()
	marginal = True  # treat insertion at the border as soft clipping
	
	for s in alnseq:  # TODO: use dict for lookup
		if s.islower() or s == "*":
			if marginal:
				cs.append("S")  # soft clipping or insertion
			else:
				cs.append("I")  # soft clipping or insertion
		elif s.isupper():
			cs.append("M")  # positional match
			marginal = False
		elif s == "-":
			cs.append("D")  # deletion
			marginal = False
	
	if len(cs) and cs._symbol[-1] == "I":
		cs._symbol[-1] = "S"  # convert trailing insertion to soft clipping
	
	return cs


def aligned_from_cigar(cigar, rawseq):
	"""Take an ungapped sequence and corresponding cigar string to create an alignment."""
	raise NotImplementedError


if __name__ == "__main__":
	from sys import stdin, stdout
	
	data = read_stockholm(stdin)
	
	for name, seq in data.items():
		stdout.write("{}\t{}\n".format(name, cigar_from_aligned_hmmalign(seq)))
