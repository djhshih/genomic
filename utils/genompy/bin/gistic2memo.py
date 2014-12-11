#!/usr/bin/env python3

import argparse, math

from genompy.cn import *


def print_peaks(peaks):
	for region in peaks:
		print(str(region))
		# display genes in region
		print(' '.join(region.genes))
		print(sum(region.samples) / float(len(region.samples)))

def create_headers(names):
	'''Create a map of header names to indices'''
	m = {}
	for i in range(len(names)):
		m[ names[i] ] = i
	return m


class RegionCoordinates:
	
	def __init__(self, peak, wpeak, region):
		self.peak = peak
		self.wpeak = wpeak
		self.region = region

	def __eq__(self, other):
		#return self.peak == other.peak and self.wpeak == other.wpeak and self.region == other.region
		return self.peak == other.peak


class CoordinatesSet:

	def __init__(self):
		self.coords = []
	
	def append(self, regionCoords):
		'''Append region coordinates, if not already in set'''
		if regionCoords not in self.coords:
			self.coords.append(regionCoords)


def main():

	pr = argparse.ArgumentParser(description='Create MEMO files from GISTIC all lesions file.')

	pr.add_argument('input', help='input GISTIC all lesions file')
	pr.add_argument('--ampfile', help='output MEMO-formated GISTIC amplified regions file', default=None)
	pr.add_argument('--delfile', help='output MEMO-formated GISTIC deleted regions file', default=None)
	pr.add_argument('--genesDb', help='knownGene sqlite database file', default='refGene.db')
	pr.add_argument('--region_attr', help='attribute name for region limits', default='Region Limits')
	pr.add_argument('--peak_attr', help='attribute name for peak limits', default='Peak Limits')
	pr.add_argument('--wpeak_attr', help='attribute name for wide peak limits', default='Wide Peak Limits')
	pr.add_argument('--delimiter', help='delimiting character for fields', default='\t')
	pr.add_argument('--item_delimiter', help='delimiting character for items within a field', default=',')

	argv = pr.parse_args()

	sets = { 'amp': CoordinatesSet(), 'del': CoordinatesSet() }

	geneDb = GeneDatabase(argv.genesDb)

	with open(argv.input, 'r') as inputf:

		header = inputf.readline()
		headers = create_headers( header.rstrip().split(argv.delimiter) )
		peak_col = headers[ argv.peak_attr ]
		wpeak_col = headers[ argv.wpeak_attr ]
		region_col = headers[ argv.region_attr ]

		for line in inputf:
			cells = line.rstrip().split(argv.delimiter)

			if 'Amplification' in cells[0]:
				state = 'amp'
			elif 'Deletion' in cells[0]:
				state = 'del'

			peak = GenomicRegion(cells[peak_col], geneDb)
			wpeak = GenomicRegion(cells[wpeak_col], geneDb)
			region = GenomicRegion(cells[region_col], geneDb)

			sets[state].append(
				RegionCoordinates(peak, wpeak, region) )

	if argv.ampfile is None:
		argv.ampfile = '{}.amp.memo'.format(argv.input)
	
	if argv.delfile is None:
		argv.delfile = '{}.del.memo'.format(argv.input)

	outputs = { 'amp': argv.ampfile, 'del': argv.delfile }

	# N.B. Last 3 columns are not used
	outheaders = [
			'index', 'chromosome',
			'region_start', 'region_end',
			'peak_start', 'peak_end',
			'enlarged_peak_start', 'enlarged_peak_end',
			'n_genes_in_region', 'genes_in_region',
			'n_genes_in_peak', 'genes_in_peak',
			'n_genes_on_chip', 'genes_on_chip', 'top 3']

	for state in outputs.keys():

		if outputs[state] is None: continue

		with open(outputs[state], 'w') as outputf:
			
			outputf.write( argv.delimiter.join(outheaders) + '\n' )

			index = 1
			for coord in sets[state].coords:

				print(str(coord.wpeak))

				wpeak_genes = coord.wpeak.genes
				region_genes = coord.region.genes

				cells = [ str(x) for x in [
					index, coord.peak.chromosome,
					coord.region.start, coord.region.end,
					coord.peak.start, coord.peak.end,
					coord.wpeak.start, coord.wpeak.end,
					len(region_genes), argv.item_delimiter.join(region_genes),
					len(wpeak_genes), argv.item_delimiter.join(wpeak_genes)
				] ]
				
				outputf.write( argv.delimiter.join(cells) + '\n' )
				index += 1


if __name__ == '__main__':
	main()

