#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.lines as lines
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

from .. import cn


def plot_sample_profile(x, y, seg_x=None, seg_y=None, subplot_spec=None, ax=None, hide_xaxis=False, yref=0, downsample=1):

	ax = plt.subplot(subplot_spec, sharex=ax) 

	# plot reference line
	ax.axhline(yref, color='#cccccc', lw=8)

	if downsample > 1:
		marker = '.'
	else:
		marker = 'o'
	
	# plot data points
	ax.plot(x[::downsample], y[::downsample], 'k', marker=marker, ls='')

	# plot segments
	if seg_x is not None and seg_y is not None:
		for i in range(len(seg_y)):

			# FIXME expose colour thresholds as parameters
			if seg_y[i] > 0.2:
				color = '#bd0026'
			elif seg_y[i] < -0.2:
				color = '#0868ac'
			else:
				color = '#666666'

			line = lines.Line2D(seg_x[i, :], [seg_y[i], seg_y[i]], color=color, lw=2)
			ax.add_line(line)

	ax.yaxis.tick_left()

	if hide_xaxis:
		for side in ('right', 'top', 'bottom'):
			ax.spines[side].set_color('none')
		ax.xaxis.set(visible=False)

	return ax


def coordinate_kbp(x, pos):
	return '{} kbp'.format(x/1e3)

def coordinate_mbp(x, pos):
	return '{} Mbp'.format(x/1e6)


def draw_xaxis(subplot_spec, ax):

	ax2 = plt.subplot(subplot_spec, sharex=ax)

	ax2.spines['right'].set_color('none')
	ax2.spines['left'].set_color('none')
	ax2.spines['top'].set_color('none')

	ax2.set(yticks=[])
	ax2.xaxis.tick_bottom()


def plot_mrna(gene_region, ax, top):
	h = 0.05

	# draw intron
	line = lines.Line2D([gene_region.start, gene_region.end], [top - h/2, top - h/2], color='#fd8d3c', lw=2, zorder=0)
	ax.add_line(line)

	# draw exons
	for exon in gene_region.exons:
		s = exon.start
		e = exon.end
		x = [s, e, e, s]
		y = [top-h, top-h, top, top]
		ax.fill(x, y, '#fd8d3c', lw=2, ec='#fd8d3c')

	# draw coding regions
	for cds in gene_region.coding_exons:
		s = cds.start
		e = cds.end
		x = [s, e, e, s]
		y = [top-h, top-h, top, top]
		ax.fill(x, y, '#f03b20', lw=2, ec='#f03b20')

	# annotate gene
	size = gene_region.end - gene_region.start + 1
	mid = gene_region.start + size/2
	ax.text(mid, -0.3, gene_region.name, horizontalalignment='right', clip_on=True, rotation=30)


def plot_mrnas(genes, subplot_spec, ax):
	ax = plt.subplot(subplot_spec, sharex=ax)

	#top = -1
	for g in genes:
		if g.strand == '+':
		#if top < 0:
			top = 0
		else:
			top = -0.1
		plot_mrna(g, ax, top)

	# remove ticks and boxes
	ax.xaxis.set(visible=False)
	ax.yaxis.set(visible=False)
	ax.set_frame_on(False)

	# set axis limits
	ax.set_ylim(bottom=-0.8, top=0)
	
	return ax


def get_gene_regions(genomicRegion):
	genes = []
	for g in genomicRegion.genes:
		genes.append(cn.GeneRegion(g, genomicRegion.geneDb))
	# sort genes by starting position
	genes.sort(key=lambda g: g.start)
	return genes


from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector, BboxConnectorPatch

def connect_bboxes(bbox1, bbox2, \
		loc1a, loc2a, loc1b, loc2b, \
		prop_lines, prop_patches=None):

	if prop_patches is None:
		prop_patches = prop_lines.copy()
		prop_patches['alpha'] = prop_patches.get('alpha', 1)*0.1
	
	c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
	c1.set_clip_on(False)
	c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
	c2.set_clip_on(False)

	bbox_patch1 = BboxPatch(bbox1, **prop_patches)
	bbox_patch2 = BboxPatch(bbox2, **prop_patches)

	p = BboxConnectorPatch(bbox1, bbox2, loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b, **prop_patches)
	p.set_clip_on(False)

	return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect(ax1, ax2, xlim, **kwargs):

	trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
	trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

	bbox = Bbox.from_extents(xlim[0], 0, xlim[1], 1)

	tbbox1 = TransformedBbox(bbox, trans1)
	tbbox2 = TransformedBbox(bbox, trans2)

	
	prop_patches = kwargs.copy()
	prop_patches['ec'] = 'none'
	prop_patches['alpha'] = 0.1

	c1, c2, bbox_patch1, bbox_patch2, p = \
			connect_bboxes(tbbox1, tbbox2, loc1a=3, loc2a=2, loc1b=4, loc2b=1, prop_lines=kwargs, prop_patches=prop_patches)
	
	ax1.add_patch(bbox_patch1)
	ax2.add_patch(bbox_patch2)
	ax2.add_patch(c1)
	ax2.add_patch(c2)
	ax2.add_patch(p)

	return c1, c2, bbox_patch1, bbox_patch2, p


def plot_locus(x, y, seg_x=None, seg_y=None, xlim=None, genes=None, yref=0, downsample=(1,1)):

	gs = gridspec.GridSpec(3, 1, height_ratios=[1, 0.75, 1])
	
	fig = plt.figure()
	fig.patch.set(facecolor='w')

	main_ax = plot_sample_profile(x, y, seg_x, seg_y, gs[0, :], yref=yref, downsample=downsample[0])
	main_ax.xaxis.set_major_formatter(ticker.FuncFormatter(coordinate_mbp))

	# zoomed-in profile
	ax = plot_sample_profile(x, y, seg_x, seg_y, gs[1, :], yref=yref, downsample=downsample[1])
	main_ax.set_xlim(left=x[0], right=x[-1])
	ax.grid(True)

	genes_ax = None
	if genes is not None:
		genes_ax = plot_mrnas(genes, gs[2, :], ax)

	if xlim is None:
		xlim = x[0], x[-1]

	win_size = xlim[1] - xlim[0] + 1
	ax.set_xlim(left=xlim[0] - win_size*0.05, right=xlim[1] + win_size*0.05)
	ax.xaxis.set_major_formatter(ticker.FuncFormatter(coordinate_mbp))

	zoom_effect(main_ax, ax, xlim, lw=0, color='green')

	return main_ax, ax, genes_ax


def main():

	geneDatabase = cn.GeneDatabase('refGene.db')
	gregion = cn.GenomicRegion('chr5:121000000-122000000', geneDatabase=geneDatabase)
	unit = 1e4

	x = np.arange(gregion.start, gregion.end+unit, unit)
	y = np.hstack( [np.random.randn(21)+1, np.random.randn(30)-1, np.random.randn(50)+2] )

	seg_x = np.matrix([ [x[0], x[20]], [x[20], x[50]], [x[50], x[-1]] ])
	seg_y = np.array([1, -1, 2])

	genes = get_gene_regions(gregion)

	xlim = (121110000, 121510000)

	main_ax, ax, genes_ax = plot_locus(x, y, seg_x, seg_y, xlim, genes)
	main_ax.set(ylabel='DNA copy-number')

	plt.show()


if __name__ == '__main__':
	main()

