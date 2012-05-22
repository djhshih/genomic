#!/usr/bin/env python3

import nose, math
from .. import cn

tolerance = 1e-6


def checkFiles(a, b):

	with open(a, 'r') as file1, open(b, 'r') as file2:

		# discard header lines
		file1.readline()
		file2.readline()

		line1 = file1.readline()
		while line1 != '':
			assert line1 == file2.readline()
			line1 = file1.readline()
		assert file2.readline() == ''


def testIO():

	s = cn.SegmentedSampleSet('test.seg')
	s.write('test.seg.out')

	checkFiles('test.seg', 'test.seg.out')


def testFind():

	s = cn.SegmentedSampleSet('test.seg')

	cnaFilter = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
	s.evaluate(cnaFilter)

	regions = [x for x in s.samples[0].find( cn.Region(chromosome='X', start=1, end=10e6) )]
	assert len(regions) == 0

	regions = [x for x in s.samples[1].find( cn.Region(chromosome='3', start=60e6, end=70e6) )]
	assert len(regions) == 2
	assert cn.Region(chromosome='3', start=26596561, end=64929919) in regions
	assert cn.Region(chromosome='3', start=64929936, end=64932005) in regions

	regions = [x for x in s.samples[2].find( cn.Region(chromosome='3', start=180366000, end=180367000) )]
	assert regions[0].state - 0.982292 < tolerance and regions[0].score == 2


