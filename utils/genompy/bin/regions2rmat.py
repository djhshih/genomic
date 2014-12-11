#!/usr/bin/env python3

import argparse, math

from genompy.cn import *


def main():

	pr = argparse.ArgumentParser(description='Create copy number data matrix (regions by samples) from a regions file.')

	pr.add_argument('input', help='input regions file')
	pr.add_argument('segfile', help='input copy number segmentation file', default=None)
	pr.add_argument('output', help='output CN file')
	pr.add_argument('--overlap', help='overlap Dice threshold', type=float, default=0.0)
	pr.add_argument('--delimiter', help='delimiting character for fields', default='\t')
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

	sampleSet = SegmentedSampleSet(argv.segfile, cnaFilter=cnaFilter)

	with open(argv.input, 'r') as inputf:

		# assume no header line

		samples = sampleSet.names

		for line in inputf:
			cells = line.rstrip().split(argv.delimiter)

			chrom = int(cells[0][ (cells[0].index('chr')+3): ])

			# add non-duplicate regions
			region = Region(chromosome=chrom, start=cells[1], end=cells[2])
			if region not in regions:
				regions.append( region )

	# populate scores matrix
	if sampleSet:

		sampleSet.prepare()

		for region in regions:
			print( str(region) )
			scores_mat.append( sampleSet.region_scores(region, argv.overlap) )


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

