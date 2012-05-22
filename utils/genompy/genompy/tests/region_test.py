#!/usr/bin/env python3

import nose
import nose.tools

from .. import cn

regionTests = [
	('chr9:12-34', cn.Region(chromosome=9, start=12, end=34)),
	('chrX:980394-1104002', cn.Region(chromosome='X', start=980394, end=1104002)),
	('chr1:5649-29485(probes 12-45)', cn.Region(chromosome=1, start=5649, end=29485)),
	('chr8:50-90(probes 11-35)|G', cn.Region(chromosome=8, start=50, end=90)),
	('chr2:294-3945|L', cn.Region(chromosome=2, start=294, end=3945)),
	('chr5', cn.Region(chromosome=5, start=0, end=0)),
]

invalidRegionTests = [
	'chr2:9-1',
	'string',
	'230',
]

def testRegion():
	for test in regionTests:
		r = cn.Region(test[0])
		print(str(r))
		assert r == test[1]
	for test in invalidRegionTests:
		nose.tools.assert_raises(test)

cytobandRegionTests = [
	('5', cn.Region(chromosome=5, start=0, end=180857866)),
	('chr9', cn.Region(chromosome=9, start=0, end=140273252)),
	('1q', cn.Region(chromosome=1, start=124300000, end=247249719)),
	('3p21', cn.Region(chromosome=3, start=43600000, end=54400000)),
	('4q31.1', cn.Region(chromosome=4, start=139500000, end=141700000)),
	('chr5q13', cn.Region(chromosome=5, start=66500000, end=76400000)),
]

def testCytobandRegion():

	cytobandTable = cn.CytobandTable('cytoBand.txt')

	for test in cytobandRegionTests:
		r = cn.CytobandRegion(test[0], cytobandTable)
		print(str(r))
		assert r == test[1]
	
