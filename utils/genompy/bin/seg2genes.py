#!/usr/bin/env python3

import argparse, os

from genompy.cn import *


pr = argparse.ArgumentParser(description='Annotate sgementation file with genes')

pr.add_argument('input', help='input regions file')
pr.add_argument('output', help='output genes file')
pr.add_argument('-d', '--delimiter', help='delimiting character', default='\t')
pr.add_argument('--item_delimiter', help='delimiting character for items within a field', default=',')
pr.add_argument('--header', help='whether to keep header', default=True)
pr.add_argument('--genesDb', help='knownGene sqlite database file', default='refGene.db')

argv = pr.parse_args()

geneDb = GeneDatabase(argv.genesDb)

regions = []
header = None
lines = []

with open(argv.input, 'r') as inputf:

	if argv.header:
		header = inputf.readline().rstrip()

	for line in inputf:

		line = line.rstrip()
		lines.append(line)

		cells = line.split(argv.delimiter)
		chrom = cells[1]
		start, end = int(cells[2]), int(cells[3])

		region = GenomicRegion(region=Region(chromosome=chrom, start=start, end=end), geneDatabase=geneDb)
		regions.append(region)


with open(argv.output, 'w') as genef:
	if header:
		genef.write( header + argv.delimiter + argv.delimiter.join(['ngenes', 'genes']) + '\n' )

	i = 0
	for r in regions:
		rs = str(r)
		print(rs)
		genes = r.genes

		genef.write( argv.delimiter.join([ lines[i], str(len(genes)), argv.item_delimiter.join(genes) ]) + '\n' )
		i += 1

