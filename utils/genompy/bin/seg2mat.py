#!/usr/bin/env python3

import argparse, math

from genompy.cn import *


def main():

	pr = argparse.ArgumentParser(description='Create gene copy number matrix from segmentation file.')

	pr.add_argument('segfile', help='input copy number segmentation file')
	pr.add_argument('output', help='output CN file')
	#pr.add_argument('--genefile', help='output genes file', default=None)
	pr.add_argument('--genesDb', help='knownGene sqlite database file', default='refGene.db')
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
	pr.add_argument('--chromosome', help='query only genes on specified chromosome', default=None)

	argv = pr.parse_args()
	if argv.cdsOnly: argv.exonsOnly = True

	regions = []
	scores_mat = []
	region_genes = {}
	sampleSet = None

	geneDb = GeneDatabase(argv.genesDb)
	cnaFilter = CNAFilter(argv.gain, argv.gainHigh, argv.gainSize,
		argv.loss, argv.lossHigh, argv.lossSize)

	sampleSet = SegmentedSampleSet(argv.segfile, cnaFilter=cnaFilter)
	samples = [x for x in sampleSet.names]

	if argv.chromosome:
		# select genes on specified chromosome
		genes = geneDb.genes_chromosome(argv.chromosome)
	else:
		# select all genes
		genes = geneDb.genome()


	# populate scores matrix

	sampleSet.prepare()

	for gene in genes:
		print(gene)
		geneRegion = GeneRegion(gene, geneDb)

		if argv.exonsOnly:
			scores = sampleSet.gene_scores(geneRegion, argv.gainCovered, argv.lossCovered, argv.cdsOnly)
		else:
			scores = sampleSet.gene_scores_relaxed(geneRegion)

		scores_mat.append( scores )

	# output scores matrix
	with open(argv.output, 'w') as outputf:

		outputf.write(argv.delimiter)
		outputf.write(argv.delimiter.join(samples) + '\n')

		for i in range(len(scores_mat)):
			outputf.write(genes[i] + argv.delimiter)	
			outputf.write(argv.delimiter.join( [ str(x) for x in scores_mat[i] ] ))
			outputf.write('\n')


if __name__ == '__main__':
	main()

