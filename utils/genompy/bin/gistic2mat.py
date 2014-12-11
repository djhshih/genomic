#!/usr/bin/env python3

import argparse, math

from genompy.cn import *

def main():

	pr = argparse.ArgumentParser(description='Create gene copy number data matrix from GISTIC all lesions file.')

	pr.add_argument('input', help='input GISTIC all lesions file')
	pr.add_argument('output', help='output CN file')
	pr.add_argument('--segfile', help='input copy number segmentation file', default=None)
	pr.add_argument('--genefile', help='output genes file', default=None)
	pr.add_argument('--genesDb', help='knownGene sqlite database file', default='refGene.db')
	pr.add_argument('--region_col', help='region column 0-index', default=2)
	pr.add_argument('--sample_col', help='first sample column 0-index', default=12)
	pr.add_argument('--delimiter', help='delimiting character for fields', default='\t')
	pr.add_argument('--item_delimiter', help='delimiting character for items within a field', default=',')
	pr.add_argument('--gain', help='copy number threshold for gain', default=0.2)
	pr.add_argument('--loss', help='copy number threshold for loss', default=-0.2)
	pr.add_argument('--gainHigh', help='copy number threshold for amplification', default=math.log(5/2, 2))
	pr.add_argument('--lossHigh', help='copy number threshold for homozygous deletion', default=math.log(0.7/2, 2))
	pr.add_argument('--gainSize', help='segment size threshold for focal gain', default=12e6)
	pr.add_argument('--lossSize', help='segment size threshold for focal loss', default=12e6)
	pr.add_argument('--exonsOnly', help='whether to only consider CNA that overlaps exons', action='store_const', const=True, default=False)
	pr.add_argument('--gainCovered', help='gene coverage threshold for gain', default=0.5)
	pr.add_argument('--lossCovered', help='gene coverage threshold = for loss', default=0.0)
	pr.add_argument('--cdsOnly', help='calculate gene coverage based on coding sequence', action='store_const', const=True, default=False)

	argv = pr.parse_args()
	if argv.cdsOnly: argv.exonsOnly = True

	regions = []
	scores_mat = []
	genes = []
	region_genes = {}
	sampleSet = None
	samples = None

	geneDb = GeneDatabase(argv.genesDb)
	cnaFilter = CNAFilter(argv.gain, argv.gainHigh, argv.gainSize,
		argv.loss, argv.lossHigh, argv.lossSize)

	if argv.segfile:
		# use segmentation file, if provided
		sampleSet = SegmentedSampleSet(argv.segfile, cnaFilter=cnaFilter)

	with open(argv.input, 'r') as inputf:

		header = inputf.readline()
		headers = header.rstrip().split(argv.delimiter)

		if sampleSet:
			samples = sampleSet.names
		else:
			# no sample set: use data from input file
			try:
				samples = headers[ argv.sample_col: ]
			except IndexError:
				raise Exception('Lesions file does not contain sample CNA status and no segmentation file specified')

		for line in inputf:
			cells = line.rstrip().split(argv.delimiter)

			state = 0
			if 'Amplification' in cells[0]:
				state = 1
			elif 'Deletion' in cells[0]:
				state = -1 

			# add non-duplicate regions
			region = AberrantRegion(cells[argv.region_col], geneDb, state)
			if region not in regions:
				if not sampleSet:
					# no sample set: use data from input file
					region.samples = [ int(x) for x in cells[ argv.sample_col: ] ]
				regions.append( region )

	# populate scores matrix
	if sampleSet:

		sampleSet.prepare()

		for region in regions:
			print( str(region) )
			if argv.exonsOnly:
				scores_map = sampleSet.scores(region, geneDb, None, argv.gainCovered, argv.lossCovered, argv.cdsOnly)
			else:
				scores_map = sampleSet.scores_relaxed(region, geneDb)

			gg = sorted(scores_map.keys())
			# add genes in region to master genes list
			genes.extend(gg)
			# append nested list from scores map to scores list-of-list
			for g in gg:
				scores_mat.append( scores_map[g] )
			region_genes[str(region)] = gg

	else:

		# use data from input file
		# N.B. the data are not completely accurate
		for region in regions:
			print( str(region) )
			# use all genes in region
			gg = region.genes
			# add genes in region to master genes list
			genes.extend(gg)
			# population matrix row
			# for each sample, entry is same for all genes in region
			for i in range(len(gg)):
				scores_mat.append(region.samples)
			region_genes[str(region)] = gg


	# output scores matrix
	with open(argv.output, 'w') as outputf:

		outputf.write(argv.delimiter)
		outputf.write(argv.delimiter.join(samples) + '\n')

		for i in range(len(scores_mat)):
			outputf.write(genes[i] + argv.delimiter)	
			outputf.write(argv.delimiter.join( [ str(x) for x in scores_mat[i] ] ))
			outputf.write('\n')

	# output genes in each lesion region
	if argv.genefile:
		with open(argv.genefile, 'w') as genef:
			for g, r in region_genes.items():
				genef.write( argv.delimiter.join([ g, str(len(r)), argv.item_delimiter.join(r) ]) + '\n' )


if __name__ == '__main__':
	main()

