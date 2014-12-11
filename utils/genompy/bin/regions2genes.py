#!/usr/bin/env python3

import argparse, os

from genompy.cn import *


pr = argparse.ArgumentParser(description='Convert regions to genes')

pr.add_argument('input', help='input regions file')
pr.add_argument('output', help='output genes file')
pr.add_argument('-d', '--delimiter', help='delimiting character', default='\t')
pr.add_argument('--item_delimiter', help='delimiting character for items within a field', default=',')
pr.add_argument('--coord_only', help='keep other fields', dest='coord_only', action='store_const', const=True, default=False)
pr.add_argument('--header', help='whether to keep header', default=False)
pr.add_argument('--genesDb', help='knownGene sqlite database file', default='refGene.db')

argv = pr.parse_args()

geneDb = GeneDatabase(argv.genesDb)

regions = []
others = {}
header = None

with open(argv.input, 'r') as inputf:

	if argv.header:
		header = inputf.readline().rstrip()

	for line in inputf:

		if argv.coord_only:
			coord = line.strip()

		else:

			idx = line.index(argv.delimiter)

			# extract coord from the first field
			coord = line[:idx]
			# copy the remaining fields into hash
			others[coord] = line[(idx+1):].rstrip()

		if coord:
			region = AberrantRegion(coord=coord, geneDatabase=geneDb)
			regions.append(region)


with open(argv.output, 'w') as genef:
	if header:
		genef.write( header + argv.delimiter + argv.delimiter.join(['n', 'genes']) + '\n' )

	for r in regions:
		rs = str(r)
		print(rs)
		genes = r.genes

		if others:
			genef.write( argv.delimiter.join([ rs, others[rs], str(len(genes)), argv.item_delimiter.join(genes) ]) + '\n' )
		else:
			genef.write( argv.delimiter.join([ rs, str(len(genes)), argv.item_delimiter.join(genes) ]) + '\n' )

