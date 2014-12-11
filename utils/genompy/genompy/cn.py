#!/usr/bin/env python3

import os, re, math, threading
import sqlite3


class GeneDatabase:

	'''Connection manager for gene annotation database'''

	def __init__(self, database='knownGene.db'):

		self.db = database

		self.conn = sqlite3.connect(':memory:')
		self.c = self.conn.cursor()

		# load database from disk into memory
		self.c.execute("attach ? as disk;", (database,) );
		if database == 'knownGene.db':
			self.c.execute("create table knownGene as select * from disk.knownGene;")
			self.c.execute("create index genePositionIdx on knownGene(name, chrom, txStart, txEnd);")
			self.c.execute("create table kgXref as select * from disk.kgXref;")
			self.c.execute("create index kgXrefIdx on kgXref(kgID)")
		elif database == 'refGene.db':
			self.c.execute("create table refGene as select * from disk.refGene;")
			self.c.execute("create index genePositionIdx on refGene(name2, chrom, txStart, txEnd);")
		else:
			raise Exception("Invalid database specified")
		self.c.execute("detach disk;")

		# commit changes
		self.conn.commit()

	def genome(self):
		'''Return list of all gene names.'''
		if self.db == 'knownGene.db':
			self.c.execute("select distinct genesymbol from kgXref where kgID in select name from knownGene")
		elif self.db == 'refGene.db':
			self.c.execute("select distinct name2 from refGene")
		else:
			raise Exception("Invalid SQL database specified")

		return [row[0] for row in self.c]

	def gene_tuple(self, gene):
		'''Return tuple(s) for gene specified by name'''
		if self.db == 'knownGene.db':
			# NB  selecting attributes from two tables require crossing, which is very inefficient in SQLite
			# TODO possible work around?
			self.c.execute("select distinct genesymbol, chrom, strand, txStart, endEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from kgXref, knownGene where kgXref.kgID == knownGene.name and kgXref.genesymbol == ?", (gene, ))
		elif self.db == 'refGene.db':
			# NB  in case of multiple refSeq transcripts, an abitrary one is returned
			# TODO resolve multiple transcript conflict
			self.c.execute("select distinct name2, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from refGene where name2 == ?", (gene,) )
		else:
			raise Exception("Invalid SQL database specified")

		return [row for row in self.c]

	def gene_tuples(self, gene):
		'''Return generator for tuple of genes specified by name.'''
		if self.db == 'knownGene.db':
			# NB  selecting attributes from two tables require crossing, which is very inefficient in SQLite
			# TODO possible work around?
			self.c.execute("select distinct genesymbol, chrom, strand, txStart, endEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from kgXref, knownGene where kgXref.kgID == knownGene.name")
		elif self.db == 'refGene.db':
			# NB  in case of multiple refSeq transcripts, an abitrary one is returned
			# TODO resolve multiple transcript conflict
			self.c.execute("select distinct name2, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds from refGene")
		else:
			raise Exception("Invalid SQL database specified")

		# return generator
		for row in self.c:
			yield row
	
	def genes(self, region):
		# execute SQL query, based on specified database
		if self.db == 'knownGene.db':
			'''
			# select genes contained within region
			self.c.execute("select distinct genesymbol from kgXref where kgID in \
				(select name from knownGene where chrom == ? and txStart >= ? and txEnd <= ?)",
				(region.chromstr, region.start, region.end))
			'''
			# select genes that overlap with region to any extent (i.e. start or end position is within region)
			self.c.execute("select distinct genesymbol from kgXref where kgID in \
				( select name from knownGene where chrom == ? and \
				  txStart <= ? and txEnd >= ? )",
				(region.chromstr, region.end, region.start))
		elif self.db == 'refGene.db':
			self.c.execute("select distinct name2 from refGene where \
				chrom == ? and txStart <= ? and txEnd >= ?",
				(region.chromstr, region.end, region.start))
		else:
			raise Exception("Invalid SQL database specified")

		return [row[0] for row in self.c]
	
	def genes_chromosome(self, chrom):
		'''Select all genes in specified chromosome'''
		# execute SQL query, based on specified database
		if self.db == 'knownGene.db':
			self.c.execute("select distinct genesymbol from kgXref where kgID in \
				( select name from knownGene where chrom == ? )",
				(region.chromstr, ))
		elif self.db == 'refGene.db':
			self.c.execute("select distinct name2 from refGene where chrom == ?",
				(region.chromstr, ))
		else:
			raise Exception("Invalid SQL database specified")

		return [row[0] for row in self.c]
	
	def __del__(self):
		self.conn.close()




class Region:

	def __init__(self, coord=None, chromosome=None, start=0, end=0, count=0, region=None):
		if coord:
			m = re.match(r'(?:chr)(?P<chrom>\d+|X|Y)(?P<coord>:(?P<start>\d+)-(?P<end>\d+))?(?P<probes>\((probes )?(?P<pstart>\d+):(?P<pend>\d+)\))?(?:\|(?P<state>.))?', coord)

			if m.group is None:
				raise Exception("Invalid region coordinate specified")

			self.chromosome = m.group('chrom')
			if m.group('coord') is not None:
				self.start = int(m.group('start'))
				self.end = int(m.group('end'))
			else:
				self.start, self.end = 0, 0
			if m.group('probes') is not None:
				self.count = int(m.group('pstart')) - int(m.group('pend')) + 1
			else:
				self.count = 0
			if m.group('state') is not None:
				state = m.group('state')
				if state == 'G':
					self.state = 1.0
				elif state == 'L':
					self.state = -1.0
				else:
					self.state = 0.0

		elif region:
			self.chromosome = region.chromosome
			self.start = region.start
			self.end = region.end
			self.count = region.count
		else:
			self.chromosome = str(chromosome)
			self.start = int(start)
			self.end = int(end)
			self.count = int(count)
		
		if self.start > self.end:
			raise Exception('Invalid region limits: {}'.format(str(self)))
	
	def __str__(self):
		return 'chr{}:{}-{}'.format(self.chromosome, self.start, self.end)

	def __eq__(self, other):
		return self.chromosome == other.chromosome and self.start == other.start and self.end == other.end

	def union(self, other):
		'''Return the union region with another region, or None if the regions do not overlap.'''
		if self.chromosome == other.chromosome and self.start <= other.end and self.end >= other.start:
			# note: do not modify self.count, since more information is necessary 
			#   to calculate the true count of the union region
			return Region(
				chromosome = self.chromosome,
				start = min(self.start, other.start),
				end = max(self.end, other.end),
				count = self.count )
		else:
			return None

	def intersect(self, other):
		'''Return the intersecting region with another region, or None if the regions do not overlap.'''
		if self.chromosome == other.chromosome and self.start <= other.end and self.end >= other.start:
			return Region(
				chromosome = self.chromosome,
				start = max(self.start, other.start),
				end = min(self.end, other.end),
				count = self.count )
		else:
			return None

	@property
	def size(self):
		return self.end - self.start + 1

	@property
	def chromstr(self):
		return 'chr{}'.format(self.chromosome)


class CytobandTable:
	
	''' Cytoband information table '''

	from collections import namedtuple
	Cytoband = namedtuple('Cytoband', ['chromosome', 'start', 'end', 'name', 'stain'])

	def __init__(self, file='cytoBand.txt'):
		cytobands = []
		# file is a UCSC cytoband file
		with open(file, 'r') as inf:
			for line in inf:
				tokens = line.strip().split('\t')
				# strip 'chr' prefix from chromosome
				tokens[0] = tokens[0][3:]
				cytobands.append( self.Cytoband(*tokens) )
		self.cytobands = cytobands	


class CytobandRegion(Region):

	''' Region specified by cytoband name '''

	def __init__(self, coord, cytobandTable):
		m = re.match(r'(?:chr)?(?P<chrom>\d+|X|Y)(?:(?P<arm>p|q)(?P<b1>\d+)?(?:\.(?P<b2>\d+)?)?)?', coord)
		# e.g., 1, 1q, chr1q, chr1q36, chr1q36.33

		if m.group is None:
			raise Exception("Invalid cytoband specified")

		chromosome = m.group('chrom')

		# construct (possibly paritial) cytoband name
		name = None
		if m.group('arm') is not None:
			name = m.group('arm')
			if m.group('b1') is not None:
				# append '.' to prevent incorrect partial matches later
				name += m.group('b1') + '.'
				if m.group('b2') is not None:
					name += m.group('b2')

		# select cytobands that match
		idx = []
		i = 0
		for cb in cytobandTable.cytobands:
			if cb.chromosome == chromosome:
				if name:
					# if name is specified, the match will be more refined
					if name == cb.name[:len(name)]:
						idx.append(i)
				else:
					idx.append(i)
			i += 1

		if len(idx) == 0:
			raise Exception("Invalid cytoband specified")

		# aggregate cytobands into one, assuming table is sorted
		# only the first index and the last matters, assuming the selected cytobands are contiguous
		start = cytobandTable.cytobands[idx[0]].start
		end = cytobandTable.cytobands[idx[-1]].end

		#super().__init__(chromosome = chromosome, start = start, end = end)
		Region.__init__(self, chromosome = chromosome, start = start, end = end)


class GeneRegion(Region):

	def __init__(self, name, geneDatabase):
		# initialize all fields required by base class
		# some of these will be overwritten below
		#super().__init__()
		Region.__init__(self)

		self.name = name
		self.geneDb = geneDatabase

		tt = geneDatabase.gene_tuple(name)
		# if there are multiple tuples, use only the first one
		# TODO allow all tuple to be considered?
		t = tt[0]

		self.chromosome, self.strand = t[1], t[2]

		# convert chromosome from 'chr?' to '?'
		if self.chromosome[0:3] == 'chr':
			self.chromosome = self.chromosome[3:]

		# set transcript start and end, use number of exons as count
		self.start, self.end, self.count = int(t[3]), int(t[4]), int(t[7])

		# create exon regions 
		exonStarts = [ int(x) for x in t[8].split(',') if x ]
		exonEnds = [ int(x) for x in t[9].split(',') if x ]

		self.exons = []

		for i in range(self.count):
			self.exons.append( Region(chromosome=self.chromosome, start=exonStarts[i], end=exonEnds[i]) )

		self.coding_exons = []

		cdsStart, cdsEnd = int(t[5]), int(t[6])
		firstCodingExon, lastCodingExon = -1, -1
		# assume exonStarts and exonEnds are in order
		# strand should not matter here
		for i in range(self.count):
			if cdsStart <= exonEnds[i]:
				# coding region has started
				firstCodingExon = i
				break

		# iterate exons in reverse order
		for i in range(self.count-1, -1, -1):
			if cdsEnd >= exonStarts[i]:
				# coding region ends in current exon
				lastCodingExon = i
				break

		# create first coding exon
		if firstCodingExon >= 0 and lastCodingExon >= 0:
			self.coding_exons.append( Region(chromosome=self.chromosome, start=cdsStart, end=exonEnds[firstCodingExon]) )
			# copy exons in between
			for i in range(firstCodingExon+1, lastCodingExon):
				self.coding_exons.append( self.exons[i] )
			# create last coding exon
			self.coding_exons.append( Region(chromosome=self.chromosome, start=exonStarts[lastCodingExon], end=cdsEnd) )

	def exon_covered(self, region):
		'''Returns the proportion of exonic regions spanned by query region'''
		overlap, total = 0, 0
		for x in self.exons:
			s = region.intersect(x)
			if s: overlap += s.size
			total += x.size
		return overlap/total

	def coding_exon_covered(self, region):
		'''Returns the proportion of coding exonic regions spanned by query region'''
		overlap, total = 0, 0
		for x in self.coding_exons:
			s = region.intersect(x)
			if s: overlap += s.size
			total += x.size
		return overlap/total


class GenomicRegion(Region):

	def __init__(self, coord=None, geneDatabase=None, region=None):
		#super().__init__(coord, region=region)
		Region.__init__(self, coord, region=region)
		self.geneDb = geneDatabase

	@property
	def genes(self):
		return self.geneDb.genes(self)


class AberrantRegion(GenomicRegion):

	def __init__(self, coord=None, geneDatabase=None, state=None, samples=None, region=None, score=0):
		#super().__init__(coord, geneDatabase, region)
		GenomicRegion.__init__(self, coord, geneDatabase, region)
		if samples:
			self.samples = samples
		else:
			self.samples = []
		self.score = int(score)

		if state is None:
			# check if self.state already exists (possible assigned by super class)
			# TODO fix bad design
			try:
				self.state
			except:
				# self.state does not already exist
				self.state = 0.0
		else:
			self.state = float(state)


	def __str__(self):

		if self.state > 0:
			state = 'G'
		elif self.state < 0:
			state = 'L'
		else:
			state = '0'

		#return '{}|{}'.format(super().__str__(), state)
		return '{}|{}'.format(Region.__str__(self), state)
	
	def delimited(self, delimiter='\t'):
		z = (str(x) for x in (self.chromosome, self.start, self.end, self.count, self.state))
		return delimiter.join(z)
	
	def evaluate(self, cnaFilter):
		self.score = cnaFilter.evaluate(self)


class CNAFilter:

	def __init__(self, gain, gainHigh, gainSize, loss, lossHigh, lossSize):
		# set thresholds
		self.gain = gain
		self.gainHigh = gainHigh
		self.gainSize = gainSize
		self.loss = loss
		self.lossHigh = lossHigh
		self.lossSize = lossSize
	
	def evaluate(self, region):
		score = 0
		if region.state > 0:
			# test region for gain
			if region.state > self.gain:
				score += 1
				if region.size < self.gainSize:
					score  +=1
					if region.state > self.gainHigh:
						score += 1
		else:
			# test region for loss
			if region.state < self.loss:
				score -= 1
				if region.size < self.lossSize:
					score -=1
					if region.state < self.lossHigh:
						score -= 1
		return score


class SegmentedSample:

	def __init__(self, name, regions=None):

		self.name = name
		if regions:
			self.regions = regions
		else:
			self.regions = []
		self.conn, self.c = None, None
	
	def append(self, region):
		'''Append region to sample regions.'''
		self.regions.append(region)

	def evaluate(self, cnaFilter):
		'''Evaluate copy number aberration scores.'''
		for region in self.regions:
			region.evaluate(cnaFilter)

	def prepare(self):
		'''Create sqlite in-memory database for sample regions, in preparation for subsequent region queries.'''
		self.conn = sqlite3.connect(':memory:')
		self.c = self.conn.cursor()	
		self.c.execute("create table regions (chrom text, start integer, end integer, score integer, state float)");
		self.c.executemany("insert into regions values(?, ?, ?, ?, ?)", self.region_tuples())
		self.c.execute("create index regionsIdx on regions(chrom, start, end);");
		self.conn.commit()

	def regions_by_chrom(self):
		chrom_regions = {}
		for region in self.regions:
			if not region.chromosome in chrom_regions:
				chrom_regions[region.chromosome] = []
			chrom_regions[region.chromosome].append(region)
		return chrom_regions

	def region_tuples(self):
		'''Return a generator for iterating region tuples.'''
		for region in self.regions:
			yield (region.chromosome, region.start, region.end, region.score, region.state)

	def find(self, query):
		'''Return sample regions that overlap with query region.'''
		if self.c is None: self.prepare()
		'''
		# select only regions that are contained or partially contained in query region
		# (i.e. region cannot completely span the query and be larger than the query)
		self.c.execute("select * from regions where \
			chrom == ? and ( (start >= ? and start <= ?) or (end >= ? and end <= ?) )",
			(query.chromosome, query.start, query.end, query.start, query.end) )
		'''
		# select regions that overlap with query region
		self.c.execute("select * from regions where \
			chrom == ? and start <= ? and end >= ?",
			(query.chromosome, query.end, query.start) )

		'''
		for row in self.c:
			yield AberrantRegion(score=row[3], state=row[4],
				region = Region(chromosome=row[0], start=row[1], end=row[2]) )
		'''
		return [ AberrantRegion(score=row[3], state=row[4],
				region = Region(chromosome=row[0], start=row[1], end=row[2]) ) for row in self.c ]
	
	def region_score(self, query, overlap_threshold=0.0):
		'''Return a score for the query region'''

		# find regions that overlap with query region
		regions = self.find(query)

		# find most extreme score
		# in case of ties, prefer losses
		score = 0
		for r in regions:
			if overlap_threshold > 0:
				overlap = query.intersect(r).size
				if overlap/query.size <= overlap_threshold:
					continue
			if abs(score) < abs(r.score):
				score = r.score
			elif abs(score) == abs(r.score):
				if r.score < score:
					score = r.score

		return score

	def gene_score(self, geneRegion, overlap_threshold_gain=0, overlap_threshold_loss=0, cds=True):
		'''Return a score for the query gene region, based on the most extreme aberrant region'''

		# find regions that overlap with query region
		regions = self.find(geneRegion)

		# use only regions that meet overlap threshold
		spanning = []
		for r in regions:
			if cds:
				coverage = geneRegion.coding_exon_covered(r)
			else:
				coverage = geneRegion.exon_covered(r)

			"""
			if geneRegion.name == "SNCAIP" and self.name == "MB-941":
				import pdb
				pdb.set_trace()
			"""

			if r.score > 0:
				if coverage > overlap_threshold_gain:
					spanning.append(r)
			elif r.score < 0:
				if coverage > overlap_threshold_loss:
					spanning.append(r)
	
		# use the most extreme value
		# in case of ties, prefer losses
		score = 0
		for r in spanning:
			if r.score != 0:
				# region is aberrant
				if abs(score) < abs(r.score):
					# save the more extreme value
					score = r.score
				elif abs(score) == abs(r.score):
					# both scores are equally extreme
					if r.score < score:
						# r.score is negative
						score = r.score

		return score

	def gene_score_relaxed(self, geneRegion):
		'''Return a score for the query gene region, based on the most extreme aberrant region'''

		# find regions that overlap with query region
		regions = self.find(geneRegion)

		# use the most extreme value
		# in case of ties, prefer losses
		score = 0
		for r in regions:
			if r.score != 0:
				# region is aberrant
				if abs(score) < abs(r.score):
					# save the more extreme value
					score = r.score
				elif abs(score) == abs(r.score):
					# both scores are equally extreme
					if r.score < score:
						# r.score is negative
						score = r.score

		return score


	def gene_score2(self, geneRegion, overlap_threshold_gain=0, overlap_threshold_loss=0, cds=True):
		'''Return a score for the query gene region, based on aberrant region with greatest overlap'''
		# TODO untested

		# find regions that overlap with query region
		regions = self.find(geneRegion)

		max_overlap_gain, max_overlap_loss = overlap_threshold_gain, overlap_threshold_loss
		region_gain, region_loss = None, None

		# find aberrant region with greatest overlap
		for r in regions:
			if r.score != 0:
				if cds:
					overlap = geneRegion.coding_exon_covered(r)
				else:
					overlap = geneRegion.exon_covered(r)

				if r.score > 0:
					if overlap > max_overlap_gain:
						max_overlap_gain = overlap
						region_gain = r
				else:
					if overlap > max_overlap_loss:
						max_overlap_loss = overlap
						region_loss = r

		# select the aberrant region with greater coverage
		overlap_prop_gain = float(max_overlap_gain) / geneRegion.size
		overlap_prop_loss = float(max_overlap_loss) / geneRegion.size

		region = None
		if region_gain:
			if region_loss:
				if region_loss >= region_gain:
					region = region_loss
			else:
				region = region_gain
		else:
			region = region_loss

		if region:
			score = region.score
		else:
			score = 0

		return score

	def scores(self, query, geneDb, geneRegions=None, overlap_threshold_gain=0, overlap_threshold_loss=0, cds=True):
		'''Return a map of genes -> scores for all genes in any unbalanced (aberrant) copy number region
		   that overlaps with the query region. If a gene has multiple scores, the most extreme score is used.
			 Only aberrant regions that cover a sufficient proportion of the gene are considered.'''

		# find all genes in query region
		if not geneRegions:
			genes = geneDb.genes(query)
			geneRegions = [ GeneRegion(gene, geneDb) for gene in genes ]

		# initialize scores to 0
		scores = {x.name : 0 for x in geneRegions}

		for geneRegion in geneRegions:
			scores[geneRegion.name] = self.gene_score(geneRegion, overlap_threshold_gain, overlap_threshold_loss, cds)

		return scores
	
	def scores_relaxed(self, query, geneDb, genes_all=None):
		'''Return a map of genes -> scores for all genes in any unbalanced (aberrant) copy number region
		   that overlaps with the query region. If a gene has multiple scores, the most extreme score is used.
			 Relaxed filter criteria: any aberrant region that overlap with gene is considered,
			 regardless of whether or how much the region overlaps with gene exons.'''

		# find all genes in query region
		if not genes_all: genes_all = geneDb.genes(query)

		# find regions that overlap with query region
		regions = self.find(query)

		# initialize scores to 0
		scores = {k : 0 for k in genes_all}

		# score genes in all regions,
		# bounded by their intersect with the query region
		for r in regions:
			if r.score != 0:
				# region is aberrant
				genes = geneDb.genes(r.intersect(query))
				for gene in genes:
					if abs(scores[gene]) < abs(r.score):
						scores[gene] = r.score
					elif abs(scores[gene]) == abs(r.score):
						if r.score < scores[gene]:
							scores[gene] = r.score
		
		return scores
	
	def __del__(self):
		if self.conn: self.conn.close()


class SegmentedSampleSet:

	delimiter = '\t'
	header = ('sample', 'chromosome', 'start', 'end', 'count', 'state')

	def __init__(self, segfile, cnaFilter=None):

		self.samples = []
		self.sample_map = {}

		with open(segfile, 'r') as f:

			# discard header
			f.readline()

			prevname = ''
			s = None
			for line in f:

				cells = line.rstrip().split(self.delimiter)
				region = AberrantRegion(state=cells[5],
					region=Region(chromosome=cells[1], start=cells[2], end=cells[3], count=cells[4]) )
				currname = cells[0]

				if currname != prevname:
					# create new sample
					s = SegmentedSample(currname)
					self.samples.append(s)
					self.sample_map[currname] = s
					prevname = currname

				# add region to the last sample in the set
				s.append(region)

		if cnaFilter: self.evaluate(cnaFilter)

	def region_scores(self, query, overlap_threshold=0):
		return [sample.region_score(query, overlap_threshold) for sample in self.samples]

	def gene_scores(self, geneRegion, overlap_threshold_gain=0, overlap_threshold_loss=0, cds=True):
		return [ sample.gene_score(geneRegion, overlap_threshold_gain, overlap_threshold_loss, cds) \
				for sample in self.samples ]

	def gene_scores_relaxed(self, geneRegion):
		return [ sample.gene_score(geneRegion) for sample in self.samples ]
	
	def scores(self, query, geneDb, genes_all=None, overlap_threshold_gain=0, overlap_threshold_loss=0, cds=True):
		# find all genes in query region (pre-compute once)
		if not genes_all: genes_all = geneDb.genes(query)

		# convert gene names to regions (pre-compute once)
		geneRegions = [ GeneRegion(gene, geneDb) for gene in genes_all ]

		# initialize the scores map-of-list
		scores = {k : [] for k in genes_all}
		
		# compute gene scores for each sample and append them to scores map
		for sample in self.samples:
			s = sample.scores(query, geneDb, geneRegions, overlap_threshold_gain, overlap_threshold_loss, cds)
			for k, v in s.items():
				scores[k].append(v)

		return scores

	def scores_relaxed(self, query, geneDb):
		# find all genes in query region (pre-compute once)
		genes_all = geneDb.genes(query)

		# initialize the scores map-of-list
		scores = {k : [] for k in genes_all}
		
		# compute gene scores for each sample and append them to scores map
		for sample in self.samples:
			s = sample.scores_relaxed(query, geneDb, genes_all)
			for k, v in s.items():
				scores[k].append(v)

		return scores


	def evaluate(self, cnaFilter):
		for sample in self.samples:
			sample.evaluate(cnaFilter)

	def write(self, segfile):
		with open(segfile, 'w') as f:
			f.write(self.delimiter.join(self.header) + '\n')
			for sample in self.samples:
				for region in sample.regions:
					f.write( sample.name + self.delimiter + region.delimited(self.delimiter) + '\n')

	def prepare(self):
		'''Create sqlite in-memory database for regions of each sample (asynchronously),
		   in preparation for subsequent region queries.'''

		tt = []
		for sample in self.samples:
			t = threading.Thread(target=sample.prepare())
			t.start()
			tt.append(t)

		for t in tt: t.join()

	def get(self, sample_name):
		'''Get segmented sample'''
		return self.sample_map[sample_name]
	
	@property
	def names(self):
		for sample in self.samples:
			yield sample.name


