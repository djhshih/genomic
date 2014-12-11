#!/usr/bin/env python2
##!/usr/bin/env python3

import os, argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py

from genompy import cn
import genompy.plot.cn as gp

matplotlib.rc('axes', lw=2)

def read_annotation(annotfile, chromfile=None, posfile=None):
	
	if os.path.exists(annotfile):
		
		print('Reading annotation from HDF5 file')
		
		with h5py.File(annotfile, 'r') as f:
			chromosomes = np.array( f['chromosome'])
			positions = np.array( f['position'] )
		
	else:
		
		with open(os.path.join(argv.root, argv.chromfile)) as chromf:
			# discard header line
			chromf.readline()
			# read chromosomes into array (of 2-character string)
			chromosomes = np.array( [ item.rstrip() for item in chromf ], dtype='S2' )
		
		with open(os.path.join(argv.root, argv.posfile)) as posf:
			# discard header line
			posf.readline()
			# read positions into array
			positions = np.array( [ int(item) for item in posf ] )
		
		# write HDF5 file for future use
		with h5py.File(argv.annotfile, 'w') as f:
			f.create_dataset('chromosome', data=chromosomes)
			f.create_dataset('position', data=positions)
		
	return (chromosomes, positions)


def read_sample(sample_path, exts):
	
	# look for data files with an appropriate extension
	sample_ext = None
	for ext in exts:
		sample_fname = '{}.{}'.format(sample_path, ext)
		if os.path.exists(sample_fname):
			sample_ext = ext
			break
			
	if sample_ext == 'hdf5':
		
		print('Reading sample profile from HDF5 file')
		
		with h5py.File(sample_fname, 'r') as f:
			profile = np.array( f[os.path.basename(sample_path)] )
		
	elif sample_ext == 'txt':
		
		with open(sample_fname, 'r') as samplef:
			sample_name = samplef.readline().strip()
			profile = np.array( [ float(item) for item in samplef ], dtype='float32' )
		
		# write HDF5 file for future use
		sample_hdf5= '{}.hdf5'.format(sample_path)
		with h5py.File(sample_hdf5, 'w') as f:
			f.create_dataset(sample_name, data=profile)
		
	else:
		raise Exception('Invalid sample file extension')
	
	return profile




pr = argparse.ArgumentParser(description='To generate genomic plot for a sample')
pr.add_argument('sample', help='name of sample')
pr.add_argument('coord', help='coordinates to plot')
pr.add_argument('coord2', help='zoomed coordinates to plot')
pr.add_argument('--seg', help='copy-number segmentation file')
pr.add_argument('--in_scale', help='original copy number scale', default='linear')
pr.add_argument('--scale', help='copy number scale', default='log')
pr.add_argument('--cytoband', help='UCSC cytoBand flat file', default='cytoBand.txt')
pr.add_argument('--geneDb', help='gene datatbase', default='refGene.db')
pr.add_argument('--annotfile', help='marker annotation file in HDF5 format', default='markers.hdf5')
#pr.add_argument('--markerfile', help='marker names file name', default='markers.txt')
pr.add_argument('--posfile', help='positions file name', default='positions.txt')
pr.add_argument('--chromfile', help='chromosomes file name', default='chromosomes.txt')
pr.add_argument('--datadir', help='data subdirectory', default='data')
pr.add_argument('--dataext', help='sample data file extension', default='hdf5,txt')
pr.add_argument('--root', help='root path to data files', default='.')
pr.add_argument('--ylim', help='y-axis limits', type=float, nargs=2, default=(-5, 5))
pr.add_argument('--downsample', help='downsampling windows', type=int, nargs=2, default=(20, 1))

argv = pr.parse_args()

scale = argv.scale
in_scale = argv.in_scale
sample_name = argv.sample

sampleSet = None
if argv.seg:
	sampleSet = cn.SegmentedSampleSet(argv.seg)

cytobandTable = cn.CytobandTable(argv.cytoband)
geneDatabase = cn.GeneDatabase(argv.geneDb)

#ylim = tuple([float(x) for x in argv.ylim])
ylim = argv.ylim
downsample = argv.downsample


annotfile, chromfile, posfile = ( os.path.join(argv.root, x) for x in (argv.annotfile, argv.chromfile, argv.posfile) )
chromosomes, positions = read_annotation(annotfile, chromfile, posfile)


sample_path = os.path.join(argv.root, argv.datadir, sample_name)
profile = read_sample(sample_path, argv.dataext.split(','))


if '-' in argv.coord:
	region = cn.Region(argv.coord)
else:
	region = cn.CytobandRegion(argv.coord, cytobandTable)


if '-' in argv.coord2:
	region2 = cn.Region(argv.coord2)
else:
	region2 = cn.CytobandRegion(argv.coord2, cytobandTable)
region2 = cn.GenomicRegion(region=region2, geneDatabase=geneDatabase)

#### Debug begin ####

'''
sample_name = 'MB-1301'
scale = 'log'
sampleSet = cn.SegmentedSampleSet('magic.seg')
cytobandTable = cn.CytobandTable('cytoBand.txt')
geneDatabase = cn.GeneDatabase('refGene.db')

chromosomes, positions = read_annotation('markers.hdf5')
profile = read_sample('data/{}'.format(sample_name), ['hdf5'])
region = cn.CytobandRegion('chr8q', cytobandTable)
region2 = cn.GenomicRegion('chr8:128700000-129500000', geneDatabase=geneDatabase)
'''

#### Debug end ####


# get sample raw copy-number estimates
idx = (chromosomes == region.chromosome.encode()) & (positions >= region.start) & (positions <= region.end)
x, y = positions[idx], profile[idx]

if scale == 'log':
	if in_scale != 'log': y = np.log2(y/2)
	yref = 0
else:
	yref = 2

seg_x = None
seg_y = None
# get sample segments
if sampleSet:
	sample = sampleSet.get(sample_name)
	sample_segs = sample.find(region)

	# create segment coordinate arrays from sample segments
	seg_x = np.zeros((len(sample_segs), 2))
	seg_y = np.zeros((len(sample_segs),))
	i = 0
	for s in sample_segs:
		seg_x[i,:] = ( (s.start, s.end) )
		seg_y[i] = s.state
		i += 1

genes = gp.get_gene_regions(region2)

# plot
main_ax, ax, genes_ax = gp.plot_locus(x, y, seg_x, seg_y, xlim=(region2.start, region2.end), genes=genes, yref=yref, downsample=downsample)
main_ax.set(ylabel='DNA copy-number')
main_ax.set(ylim=ylim)
ax.set(ylabel='DNA copy-number')
ax.set(ylim=ylim)

plt.show()

