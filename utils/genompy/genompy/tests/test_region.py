#!/usr/bin/env python3

import os, unittest, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import cn

test_dir = os.path.dirname(os.path.abspath(__file__))


class RegionTest(unittest.TestCase):

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

    def testRegion(self):
        old = os.getcwd()
        os.chdir(test_dir)
        try:
            for coord, expected in self.regionTests:
                r = cn.Region(coord)
                self.assertEqual(r, expected)
            for coord in self.invalidRegionTests:
                self.assertRaises(Exception, cn.Region, coord)
        finally:
            os.chdir(old)


class CytobandRegionTest(unittest.TestCase):

    cytobandRegionTests = [
        ('5', cn.Region(chromosome=5, start=0, end=180857866)),
        ('chr9', cn.Region(chromosome=9, start=0, end=140273252)),
        ('1q', cn.Region(chromosome=1, start=124300000, end=247249719)),
        ('3p21', cn.Region(chromosome=3, start=43600000, end=54400000)),
        ('4q31.1', cn.Region(chromosome=4, start=139500000, end=141700000)),
        ('chr5q13', cn.Region(chromosome=5, start=66500000, end=76400000)),
    ]

    def testCytobandRegion(self):
        old = os.getcwd()
        os.chdir(test_dir)
        try:
            ct = cn.CytobandTable('cytoBand.txt')
            for coord, expected in self.cytobandRegionTests:
                r = cn.CytobandRegion(coord, ct)
                self.assertEqual(r, expected)
        finally:
            os.chdir(old)


if __name__ == '__main__':
    unittest.main()
