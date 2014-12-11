#!/usr/bin/env python3

import os, argparse

from cn.cn import *

pr = argparse.ArgumentParser('To annote genes table with coordinate information')
pr.add_argument('input', help='input genes table file')
pr.add_argument('output', help='output genes table file')
pr.add_argument('--delimiter', help='delimiting character', default='\t')
pr.add_argument('--genesDb', help='genes annotation sqlite data base', default='refGene.db')

argv = pr.parse_args()

geneDb = GeneDatabase(argv.genesDb)

header = None
lines = []
geneRegions = []

with open(argv.input, 'r') as inputf:

	header = inputf.readline().rstrip()

	for line in inputf:

		line = line.strip()
		lines.append(line)

		cells = line.split(argv.delimiter)
		gene = cells[0]

		geneRegions.append( GeneRegion(gene, geneDb) )


with open(argv.output, 'w') as outputf:

	outputf.write('{}{}{}\n'.format( header, argv.delimiter, argv.delimiter.join(['chromosome', 'start', 'end']) ))

	for i in range(len(geneRegions)):

		geneRegion = geneRegions[i]

		coords = [ str(x) for x in [geneRegion.chromosome, geneRegion.start, geneRegion.end] ]
		outputf.write( '{}{}{}\n'.format(lines[i], argv.delimiter, argv.delimiter.join(coords)) )

