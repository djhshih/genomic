#!/usr/bin/env python3

import random

from . import cn
from . import gp

def overlap_genes_in_regions(regions, geneDb, geneSets, overlaps, genes=None):
	'''Track overlap of genes in regions with genes in each gene set, in place'''
	# genes and overlaps will be modified in place

	# store all genes from all regions
	if genes is None:
		genes = set()

	for region in regions:
		region_genes = geneDb.genes(region)
		for gene in region_genes:
			genes.add(gene)

	overlap_genes_in_set(genes, geneSets, overlaps)

def overlap_genes_in_set(genes, geneSets, overlaps):
	'''Track overlap of genes with each gene set, in place'''
	# overlaps will be modified in place
	
	# determine overlap of genes with gene sets

	for name, gs in geneSets.sets.items():
		gene_set = gs.genes

		# count overlap
		overlap = 0
		for gene in gene_set:
			if gene in genes:
				overlap += 1

		overlaps[name].append(overlap)

def randomized_genes(genes, universe):
	'''Return set of re-sampled genes from universe'''
	new_genes = set()
	for i in range(len(genes)):
		new_genes.add( random.choice(universe) )
	return new_genes

def randomize_regions(regions, chrlens, rand_chrom=True):
	'''Randomize regions in place'''
	for region in regions:
		randomize_region(region, chrlens, rand_chrom)

def randomize_region(region, chrlens, rand_chrom=True):
	'''Randomize region in place'''
	if rand_chrom:
		# randomly choose a new chromosome
		region.chromosome = random.choice([x for x in chrlens.keys()])
	
	# randomize start position
	size = region.size
	max_end = chrlens[region.chromosome] - size + 1
	region.start = random.randint(0, max_end)
	region.end = region.start + size - 1

