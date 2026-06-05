#!/usr/bin/env python3

import os, math, unittest, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import cn

tolerance = 1e-6
test_dir = os.path.dirname(os.path.abspath(__file__))


class CNTest(unittest.TestCase):

    def checkFiles(self, a, b):
        with open(a, 'r') as file1, open(b, 'r') as file2:
            file1.readline()
            file2.readline()
            line1 = file1.readline()
            while line1 != '':
                self.assertEqual(line1, file2.readline())
                line1 = file1.readline()
            self.assertEqual(file2.readline(), '')

    def testIO(self):
        old = os.getcwd()
        os.chdir(test_dir)
        try:
            s = cn.SegmentedSampleSet('test.seg')
            s.write('test.seg.out')
            self.checkFiles('test.seg', 'test.seg.out')
        finally:
            os.chdir(old)

    def testFind(self):
        old = os.getcwd()
        os.chdir(test_dir)
        try:
            s = cn.SegmentedSampleSet('test.seg')
            f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6,
                             -0.1, math.log(0.7/2, 2), 12e6)
            s.evaluate(f)

            regions = s.samples[0].find(cn.Region(chromosome='X', start=1, end=10e6))
            self.assertEqual(len(regions), 0)

            regions = s.samples[1].find(cn.Region(chromosome='3', start=60e6, end=70e6))
            self.assertEqual(len(regions), 2)
            self.assertIn(cn.Region(chromosome='3', start=26596561, end=64929919), regions)
            self.assertIn(cn.Region(chromosome='3', start=64929936, end=64932005), regions)

            regions = s.samples[2].find(cn.Region(chromosome='3', start=180366000, end=180367000))
            self.assertTrue(regions[0].state - 0.982292 < tolerance)
            self.assertEqual(regions[0].score, 2)
        finally:
            os.chdir(old)


if __name__ == '__main__':
    unittest.main()
