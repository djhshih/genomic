#!/usr/bin/env python3

import os

class GeneSet:
	
	def __init__(self, name='', description='', genes=None, txt_file=None):
		self.name = name
		self.description = description
		if genes is None:
			self.genes = set()
		else:
			self.genes = set(genes)

		if txt_file:
			self.read(txt_file)

	def read(self, txt_file, derive_name=False):
		'''Read in genes from txt file, containing one gene per row'''
		if derive_name:
			# set name to file basename without extension
			s = os.path.basename(txt_file)
			idx = s.rfind('.')
			if idx != -1:
				s = s[:idx]
			self.name = s

		# read in genes, one per row
		with open(txt_file, 'r') as f:
			for line in f:
				self.genes.add(line.strip())

	def __str__(self):
		return self.name


class GeneSets:

	delimiter = '\t'

	def __init__(self, gmt_file=None):
		self.sets = {}

		if gmt_file:
			self.read(gmt_file)
	
	def read(self, gmt_file):
		'''Read gene sets from gmt file'''
		with open(gmt_file, 'r') as f:
			for line in f:
				cells = line.strip().split(self.delimiter)
				# first cell is name, second is description, rest are genes
				gs = GeneSet(name=cells[0], description=cells[1], genes=cells[2:])
				self.sets[gs.name] = gs
	
	def write(self, gmt_file):
		'''Write gene sets to gmt file'''
		with open(gmt_file, 'w') as f:
			for k, s in self.sets.items():
				t = [s.name, s.description]
				t.extend(s.genes)
				line = self.delimiter.join(t)
				f.write(line + '\n')
	
	def add(self, gene_set):
		'''Add a gene set to sets'''
		self.sets[gene_set.name] = gene_set

