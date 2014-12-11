#!/usr/bin/env python3

import argparse, math

from genompy.cn import *


def main():

	pr = argparse.ArgumentParser(description='Create copy number data matrix (regions by samples) from a GISTIC all lesions file.')

	pr.add_argument('input', help='input GISTIC all lesions file')
	pr.add_argument('output', help='output CN file')
	pr.add_argument('--overlap', help='query overlap threshold', type=float, default=0.0)
	pr.add_argument('--segfile', help='input copy number segmentation file', default=None)
	pr.add_argument('--delimiter', help='delimiting character for fields', default='\t')
	pr.add_argument('--region_col', help='region column 0-index', default=2)
	pr.add_argument('--sample_col', help='first sample column 0-index', default=12)
	pr.add_argument('--gain', help='copy number threshold for gain', default=0.2)
	pr.add_argument('--loss', help='copy number threshold for loss', default=-0.2)
	pr.add_argument('--gainHigh', help='copy number threshold for amplification', default=math.log(5/2, 2))
	pr.add_argument('--lossHigh', help='copy number threshold for homozygous deletion', default=math.log(0.7/2, 2))
	pr.add_argument('--gainSize', help='segment size threshold for focal gain', default=12e6)
	pr.add_argument('--lossSize', help='segment size threshold for focal loss', default=12e6)

	argv = pr.parse_args()

	regions = []
	scores_mat = []
	genes = []
	region_genes = {}
	sampleSet = None
	samples = None

	cnaFilter = CNAFilter(argv.gain, argv.gainHigh, argv.gainSize,
		argv.loss, argv.lossHigh, argv.lossSize)

	if argv.segfile:
		sampleSet = SegmentedSampleSet(argv.segfile, cnaFilter=cnaFilter)

	with open(argv.input, 'r') as inputf:

		# assume no header line
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
			region = AberrantRegion(cells[argv.region_col], state=state)
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
			scores_mat.append( sampleSet.region_scores(region, argv.overlap) )
	
	else:
		for region in regions:
			print( str(region) )
			scores_mat.append(region.samples)


	# output scores matrix
	with open(argv.output, 'w') as outputf:

		outputf.write(argv.delimiter)
		outputf.write(argv.delimiter.join(samples) + '\n')

		for i in range(len(scores_mat)):
			outputf.write(str(regions[i]) + argv.delimiter)	
			outputf.write(argv.delimiter.join( [ str(x) for x in scores_mat[i] ] ))
			outputf.write('\n')


if __name__ == '__main__':
	main()

