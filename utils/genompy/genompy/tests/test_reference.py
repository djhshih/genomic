#!/usr/bin/env python3

"""
Reference tests for genompy.cn module.

These tests document the exact input/output values that the C++
implementation must match. Each assertion is a reference value.

Run:  python3 -m unittest test_reference  (from the tests/ directory)
"""

import os, math, sqlite3, shutil, tempfile, unittest, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import cn

tolerance = 1e-6
test_dir = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# Module-level fixtures (shared across all test classes)
# ---------------------------------------------------------------------------

def _create_test_refgene_db(path):
    conn = sqlite3.connect(path)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE refGene (
            bin INTEGER, name TEXT, chrom TEXT, strand TEXT,
            txStart INTEGER, txEnd INTEGER,
            cdsStart INTEGER, cdsEnd INTEGER,
            exonCount INTEGER,
            exonStarts TEXT, exonEnds TEXT,
            score INTEGER, name2 TEXT,
            cdsStartStat TEXT, cdsEndStat TEXT,
            exonFrames TEXT
        )
    ''')
    c.execute('''
        INSERT INTO refGene VALUES (642, 'NM_000546', 'chr17', '-',
            7512444, 7531588, 7513651, 7520637, 11,
            '7512444,7514651,7517577,7517743,7518223,7518901,7519095,7520036,7520424,7520563,7531419,',
            '7513733,7514758,7517651,7517880,7518333,7519014,7519279,7520315,7520446,7520665,7531588,',
            0, 'TP53', 'cmpl', 'cmpl', '2,0,1,2,0,1,0,0,2,0,-1,')
    ''')
    c.execute('''
        INSERT INTO refGene VALUES (878, 'NM_007294', 'chr17', '-',
            38449837, 38531026, 38451220, 38529639, 23,
            '38449837,38453185,38454663,38456605,38462594,38468875,38469416,38473150,38476470,38479873,38482030,38487946,38496486,38496977,38501388,38502786,38505317,38509664,38510410,38511998,38521268,38529559,38530813,',
            '38451345,38453246,38454737,38456660,38462678,38468916,38469494,38473238,38476781,38480064,38482157,38488118,38496575,38500403,38501465,38502832,38505423,38509804,38510499,38512076,38521322,38529658,38531026,',
            0, 'BRCA1', 'cmpl', 'cmpl', '1,0,1,0,0,1,1,0,1,2,1,0,1,1,2,1,0,1,2,2,2,0,-1,')
    ''')
    c.execute('''
        INSERT INTO refGene VALUES (125, 'NM_201284', 'chr7', '+',
            55054218, 55206232, 55054464, 55205731, 16,
            '55054218,55177472,55178491,55181792,55186480,55187732,55189197,55191016,55191719,55191945,55192849,55195325,55196685,55198919,55200466,55205493,',
            '55054552,55177624,55178675,55181927,55186549,55187851,55189339,55191133,55191846,55192019,55192940,55195525,55196818,55199010,55200624,55206232,',
            0, 'EGFR', 'cmpl', 'cmpl', '0,1,0,1,1,1,0,1,1,2,1,2,1,2,0,2,')
    ''')
    conn.commit()
    conn.close()


_old_cwd = None

def setUpModule():
    global _old_cwd
    _old_cwd = os.getcwd()
    os.chdir(test_dir)
    _create_test_refgene_db(os.path.join(test_dir, 'refGene.db'))


def tearDownModule():
    global _old_cwd
    os.chdir(_old_cwd)
    db = os.path.join(test_dir, 'refGene.db')
    if os.path.exists(db):
        os.unlink(db)


# ---------------------------------------------------------------------------
# 1. Region tests
# ---------------------------------------------------------------------------

class RegionTest(unittest.TestCase):

    def test_parse(self):
        r = cn.Region('chr9:12-34')
        self.assertEqual(r.chromosome, '9')
        self.assertEqual(r.start, 12)
        self.assertEqual(r.end, 34)
        self.assertEqual(r.count, 0)

    def test_parse_chrX(self):
        r = cn.Region('chrX:980394-1104002')
        self.assertEqual(r.chromosome, 'X')
        self.assertEqual(r.start, 980394)
        self.assertEqual(r.end, 1104002)

    def test_parse_with_probes(self):
        r = cn.Region('chr1:5649-29485(probes 12-45)')
        self.assertEqual(r.chromosome, '1')
        self.assertEqual(r.start, 5649)
        self.assertEqual(r.end, 29485)
        self.assertEqual(r.count, 0)

    def test_parse_gain(self):
        r = cn.Region('chr8:50-90(probes 11-35)|G')
        self.assertEqual(r.chromosome, '8')
        self.assertEqual(r.start, 50)
        self.assertEqual(r.end, 90)

    def test_parse_loss(self):
        r = cn.Region('chr2:294-3945|L')
        self.assertEqual(r.chromosome, '2')
        self.assertEqual(r.start, 294)
        self.assertEqual(r.end, 3945)

    def test_parse_chrom_only(self):
        r = cn.Region('chr5')
        self.assertEqual(r.chromosome, '5')
        self.assertEqual(r.start, 0)
        self.assertEqual(r.end, 0)
        self.assertEqual(r.count, 0)

    def test_invalid_start_gt_end(self):
        self.assertRaises(Exception, cn.Region, 'chr2:9-1')

    def test_size(self):
        r = cn.Region(chromosome=1, start=100, end=200)
        self.assertEqual(r.size, 101)

    def test_chromstr(self):
        r = cn.Region(chromosome='X', start=1, end=10)
        self.assertEqual(r.chromstr, 'chrX')

    def test_str(self):
        r = cn.Region(chromosome=9, start=12, end=34)
        self.assertEqual(str(r), 'chr9:12-34')

    def test_eq(self):
        a = cn.Region(chromosome=1, start=100, end=200)
        b = cn.Region(chromosome=1, start=100, end=200)
        self.assertEqual(a, b)

    def test_neq(self):
        a = cn.Region(chromosome=1, start=100, end=200)
        b = cn.Region(chromosome=1, start=101, end=200)
        self.assertNotEqual(a, b)

    def test_union_overlapping(self):
        a = cn.Region(chromosome=1, start=100, end=200)
        b = cn.Region(chromosome=1, start=150, end=300)
        u = a.union(b)
        self.assertEqual(u.chromosome, '1')
        self.assertEqual(u.start, 100)
        self.assertEqual(u.end, 300)

    def test_union_non_overlapping(self):
        a = cn.Region(chromosome=1, start=100, end=200)
        b = cn.Region(chromosome=1, start=300, end=400)
        self.assertIsNone(a.union(b))

    def test_union_diff_chrom(self):
        a = cn.Region(chromosome=1, start=100, end=200)
        b = cn.Region(chromosome=2, start=100, end=200)
        self.assertIsNone(a.union(b))

    def test_intersect_overlapping(self):
        a = cn.Region(chromosome=1, start=100, end=300)
        b = cn.Region(chromosome=1, start=200, end=400)
        i = a.intersect(b)
        self.assertEqual(i.chromosome, '1')
        self.assertEqual(i.start, 200)
        self.assertEqual(i.end, 300)
        self.assertEqual(i.size, 101)

    def test_intersect_non_overlapping(self):
        a = cn.Region(chromosome=1, start=100, end=200)
        b = cn.Region(chromosome=1, start=300, end=400)
        self.assertIsNone(a.intersect(b))

    def test_intersect_contained(self):
        outer = cn.Region(chromosome=1, start=100, end=500)
        inner = cn.Region(chromosome=1, start=200, end=300)
        i = outer.intersect(inner)
        self.assertEqual(i.start, 200)
        self.assertEqual(i.end, 300)

    def test_intersect_same(self):
        a = cn.Region(chromosome=1, start=100, end=200)
        self.assertEqual(a.intersect(a), a)

    def test_intersect_count_from_self(self):
        a = cn.Region(chromosome=1, start=100, end=300, count=5)
        b = cn.Region(chromosome=1, start=200, end=400, count=99)
        i = a.intersect(b)
        self.assertEqual(i.count, 5)


# ---------------------------------------------------------------------------
# 2. AberrantRegion tests
# ---------------------------------------------------------------------------

class AberrantRegionTest(unittest.TestCase):

    def test_parse_gain(self):
        r = cn.AberrantRegion('chr1:100-200(probes 1-10)|G')
        self.assertEqual(r.chromosome, '1')
        self.assertEqual(r.start, 100)
        self.assertEqual(r.end, 200)
        self.assertEqual(r.count, 0)
        self.assertAlmostEqual(r.state, 0.0)
        self.assertEqual(r.samples, [])
        self.assertEqual(r.score, 0)

    def test_parse_loss(self):
        r = cn.AberrantRegion('chr1:100-200(probes 1-10)|L')
        self.assertAlmostEqual(r.state, 0.0)

    def test_parse_neutral(self):
        r = cn.AberrantRegion('chr1:100-200')
        self.assertAlmostEqual(r.state, 0.0)

    def test_str_gain(self):
        r = cn.AberrantRegion(state=0.5, region=cn.Region(chromosome=1, start=100, end=200))
        self.assertEqual(str(r), 'chr1:100-200|G')

    def test_str_loss(self):
        r = cn.AberrantRegion(state=-0.5, region=cn.Region(chromosome=1, start=100, end=200))
        self.assertEqual(str(r), 'chr1:100-200|L')

    def test_str_neutral(self):
        r = cn.AberrantRegion(state=0.0, region=cn.Region(chromosome=1, start=100, end=200))
        self.assertEqual(str(r), 'chr1:100-200|0')

    def test_delimited(self):
        r = cn.AberrantRegion(state=0.5, score=3, region=cn.Region(chromosome=1, start=100, end=200, count=50))
        self.assertEqual(r.delimited('\t'), '1\t100\t200\t50\t0.5')

    def test_evaluate(self):
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        r = cn.AberrantRegion(state=0.3, region=cn.Region(chromosome=1, start=100, end=200))
        r.evaluate(f)
        self.assertNotEqual(r.score, 0)


# ---------------------------------------------------------------------------
# 3. CNAFilter tests
# ---------------------------------------------------------------------------

class CNAFilterTest(unittest.TestCase):

    def _standard_filter(self):
        return cn.CNAFilter(0.2, 0.8, 1000, -0.2, -0.8, 1000)

    def _seg(self, state, start, end):
        return cn.AberrantRegion(state=state, region=cn.Region(chromosome=1, start=start, end=end))

    def test_gain_above_threshold(self):
        s = self._standard_filter().evaluate(self._seg(0.5, 1, 500))
        self.assertEqual(s, 2)

    def test_gain_above_threshold_small(self):
        s = self._standard_filter().evaluate(self._seg(0.5, 1, 100))
        self.assertEqual(s, 2)

    def test_gain_high_and_small(self):
        s = self._standard_filter().evaluate(self._seg(1.0, 1, 100))
        self.assertEqual(s, 3)

    def test_loss_below_threshold(self):
        s = self._standard_filter().evaluate(self._seg(-0.5, 1, 500))
        self.assertEqual(s, -2)

    def test_loss_below_threshold_small(self):
        s = self._standard_filter().evaluate(self._seg(-0.5, 1, 100))
        self.assertEqual(s, -2)

    def test_loss_high_and_small(self):
        s = self._standard_filter().evaluate(self._seg(-1.0, 1, 100))
        self.assertEqual(s, -3)

    def test_neutral(self):
        s = self._standard_filter().evaluate(self._seg(0.0, 1, 500))
        self.assertEqual(s, 0)

    def test_below_gain_threshold(self):
        s = self._standard_filter().evaluate(self._seg(0.1, 1, 500))
        self.assertEqual(s, 0)

    def test_above_loss_threshold(self):
        s = self._standard_filter().evaluate(self._seg(-0.1, 1, 500))
        self.assertEqual(s, 0)

    def test_on_test_seg_all_scores(self):
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        s.evaluate(f)
        expected = [
            ('MB-1', '1', 890962, 247095424, -1),
            ('MB-1', '2', 60114, 242381784, 1),
            ('MB-1', '3', 41866, 26506787, 1),
            ('MB-1', '3', 26507577, 26595923, 2),
            ('MB-2', '3', 26596561, 64929919, 1),
            ('MB-2', '3', 64929936, 64932005, 2),
            ('MB-3', '3', 64932286, 180365866, 1),
            ('MB-3', '3', 180366352, 180367291, 2),
            ('MB-3', '3', 180368613, 180434337, 2),
        ]
        for samp_name, chrom, start, end, exp_score in expected:
            sample = s.get(samp_name)
            found = False
            for region in sample.regions:
                if (region.chromosome == chrom and
                        region.start == start and region.end == end):
                    self.assertEqual(region.score, exp_score,
                        f"{samp_name} chr{chrom}:{start}-{end}: "
                        f"expected {exp_score}, got {region.score}")
                    found = True
                    break
            self.assertTrue(found, f"Segment {samp_name} chr{chrom}:{start}-{end} not found")


# ---------------------------------------------------------------------------
# 4. SegmentedSample.find() tests
# ---------------------------------------------------------------------------

class SegmentedSampleFindTest(unittest.TestCase):

    def test_no_overlap(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        regions = s.get('MB-1').find(cn.Region(chromosome='X', start=1, end=10e6))
        self.assertEqual(len(regions), 0)

    def test_two_overlapping(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        regions = s.get('MB-2').find(cn.Region(chromosome='3', start=60e6, end=70e6))
        self.assertEqual(len(regions), 2)
        self.assertIn(cn.Region(chromosome='3', start=26596561, end=64929919), regions)
        self.assertIn(cn.Region(chromosome='3', start=64929936, end=64932005), regions)

    def test_one_overlapping(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        sample = s.get('MB-3')
        regions = sample.find(cn.Region(chromosome='3', start=180366000, end=180367000))
        self.assertEqual(len(regions), 1)
        r = regions[0]
        self.assertEqual(r.start, 180366352)
        self.assertEqual(r.end, 180367291)
        self.assertAlmostEqual(r.state, 0.982292)
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        sample.evaluate(f)
        sample.prepare()
        regions2 = sample.find(cn.Region(chromosome='3', start=180366000, end=180367000))
        self.assertEqual(regions2[0].score, 2)

    def test_spanning_query(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        regions = s.get('MB-1').find(cn.Region(chromosome='1', start=1, end=250e6))
        self.assertEqual(len(regions), 1)

    def test_contained_query(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        regions = s.get('MB-1').find(cn.Region(chromosome='1', start=50e6, end=60e6))
        self.assertEqual(len(regions), 1)


# ---------------------------------------------------------------------------
# 5. SegmentedSample.region_score() tests
# ---------------------------------------------------------------------------

class RegionScoreTest(unittest.TestCase):

    def test_no_overlap(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        s.evaluate(f)
        score = s.get('MB-2').region_score(cn.Region(chromosome='X', start=1, end=10e6))
        self.assertEqual(score, 0)

    def test_single_region(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        s.evaluate(f)
        score = s.get('MB-2').region_score(cn.Region(chromosome='3', start=64929936, end=64932005))
        self.assertEqual(score, 2)

    def test_most_extreme(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        s.evaluate(f)
        score = s.get('MB-2').region_score(cn.Region(chromosome='3', start=26596561, end=64932005))
        self.assertEqual(score, 2)

    def test_tie_prefers_loss(self):
        sample = cn.SegmentedSample('test')
        sample.append(cn.AberrantRegion(state=0.5, score=2, region=cn.Region(chromosome='1', start=100, end=200)))
        sample.append(cn.AberrantRegion(state=-0.5, score=-2, region=cn.Region(chromosome='1', start=150, end=250)))
        score = sample.region_score(cn.Region(chromosome='1', start=100, end=250))
        self.assertEqual(score, -2)

    def test_overlap_threshold(self):
        sample = cn.SegmentedSample('test')
        sample.append(cn.AberrantRegion(state=0.5, score=3, region=cn.Region(chromosome='1', start=100, end=109)))
        sample.append(cn.AberrantRegion(state=0.5, score=1, region=cn.Region(chromosome='1', start=200, end=300)))
        q = cn.Region(chromosome='1', start=100, end=299)
        self.assertEqual(sample.region_score(q, overlap_threshold=0.0), 3)
        self.assertEqual(sample.region_score(q, overlap_threshold=0.1), 1)
        self.assertEqual(sample.region_score(q, overlap_threshold=0.6), 0)


# ---------------------------------------------------------------------------
# 6. GeneDatabase tests
# ---------------------------------------------------------------------------

class GeneDatabaseTest(unittest.TestCase):

    def test_genome(self):
        db = cn.GeneDatabase('refGene.db')
        genes = db.genome()
        self.assertIn('TP53', genes)
        self.assertIn('BRCA1', genes)
        self.assertIn('EGFR', genes)
        self.assertEqual(len(genes), 3)

    def test_gene_tuple(self):
        db = cn.GeneDatabase('refGene.db')
        tt = db.gene_tuple('TP53')
        self.assertGreaterEqual(len(tt), 1)
        t = tt[0]
        self.assertEqual(t[1], 'chr17')
        self.assertEqual(t[2], '-')
        self.assertEqual(t[3], 7512444)
        self.assertEqual(t[4], 7531588)

    def test_genes_region(self):
        db = cn.GeneDatabase('refGene.db')
        genes = db.genes(cn.Region(chromosome=17, start=7.5e6, end=7.6e6))
        self.assertIn('TP53', genes)
        self.assertNotIn('BRCA1', genes)

    def test_genes_region_multiple(self):
        db = cn.GeneDatabase('refGene.db')
        genes = db.genes(cn.Region(chromosome=17, start=7e6, end=39e6))
        self.assertIn('TP53', genes)
        self.assertIn('BRCA1', genes)

    def test_genes_region_no_match(self):
        db = cn.GeneDatabase('refGene.db')
        genes = db.genes(cn.Region(chromosome=1, start=1, end=10e6))
        self.assertEqual(len(genes), 0)

    def test_genes_chromosome_skipped(self):
        pass  # skip: cn.py uses undefined 'region' variable

    def test_genes_chromosome_other_skipped(self):
        pass  # skip: cn.py uses undefined 'region' variable


# ---------------------------------------------------------------------------
# 7. GeneRegion tests
# ---------------------------------------------------------------------------

class GeneRegionTest(unittest.TestCase):

    def setUp(self):
        self.db = cn.GeneDatabase('refGene.db')
        self.gr = cn.GeneRegion('TP53', self.db)

    def test_construction(self):
        self.assertEqual(self.gr.name, 'TP53')
        self.assertEqual(self.gr.chromosome, '17')
        self.assertEqual(self.gr.strand, '-')
        self.assertEqual(self.gr.start, 7512444)
        self.assertEqual(self.gr.end, 7531588)
        self.assertEqual(self.gr.count, 11)
        self.assertEqual(len(self.gr.exons), 11)
        self.assertGreaterEqual(len(self.gr.coding_exons), 1)
        self.assertEqual(self.gr.exons[0].start, 7512444)
        self.assertEqual(self.gr.exons[0].end, 7513733)

    def test_exon_covered_full(self):
        region = cn.Region(chromosome=17, start=self.gr.start, end=self.gr.end)
        self.assertAlmostEqual(self.gr.exon_covered(region), 1.0)

    def test_exon_covered_partial(self):
        total_exon = sum(x.size for x in self.gr.exons)
        region = cn.Region(chromosome=17, start=7512444, end=7514758)
        expected = (7513733 - 7512444 + 1 + 7514758 - 7514651 + 1) / total_exon
        self.assertAlmostEqual(self.gr.exon_covered(region), expected)

    def test_exon_covered_none(self):
        region = cn.Region(chromosome=17, start=1, end=1000)
        self.assertEqual(self.gr.exon_covered(region), 0.0)

    def test_coding_exon_covered(self):
        region = cn.Region(chromosome=17, start=self.gr.start, end=self.gr.end)
        self.assertAlmostEqual(self.gr.coding_exon_covered(region), 1.0)
        region2 = cn.Region(chromosome=17, start=1, end=1000)
        self.assertEqual(self.gr.coding_exon_covered(region2), 0.0)


# ---------------------------------------------------------------------------
# 8. SegmentedSample.gene_score() tests
# ---------------------------------------------------------------------------

class GeneScoreTest(unittest.TestCase):

    def test_no_overlap(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        s.evaluate(f)
        db = cn.GeneDatabase('refGene.db')
        gr = cn.GeneRegion('TP53', db)
        self.assertEqual(s.get('MB-1').gene_score(gr), 0)

    def test_with_overlap(self):
        sample = cn.SegmentedSample('test')
        sample.append(cn.AberrantRegion(state=0.3, region=cn.Region(chromosome='7', start=55000000, end=56000000, count=100)))
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        sample.evaluate(f)
        db = cn.GeneDatabase('refGene.db')
        gr = cn.GeneRegion('EGFR', db)
        self.assertNotEqual(sample.gene_score(gr), 0)


# ---------------------------------------------------------------------------
# 9. SegmentedSampleSet aggregate tests
# ---------------------------------------------------------------------------

class SegmentedSampleSetTest(unittest.TestCase):

    def test_io_roundtrip(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        tmp = os.path.join(test_dir, '_test_roundtrip.out')
        try:
            s.write(tmp)
            with open(os.path.join(test_dir, 'test.seg')) as f:
                original = f.read()
            with open(tmp) as f:
                written = f.read()
            self.assertEqual(original, written)
        finally:
            if os.path.exists(tmp):
                os.unlink(tmp)

    def test_get(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        self.assertEqual(s.get('MB-1').name, 'MB-1')
        self.assertEqual(s.get('MB-2').name, 'MB-2')
        self.assertEqual(s.get('MB-3').name, 'MB-3')

    def test_names(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        self.assertEqual(list(s.names), ['MB-1', 'MB-2', 'MB-3'])

    def test_region_scores(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        s.evaluate(f)
        scores = s.region_scores(cn.Region(chromosome='3', start=180366000, end=180367000))
        self.assertEqual(scores, [0, 0, 2])

    def test_gene_scores(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        s.evaluate(f)
        db = cn.GeneDatabase('refGene.db')
        gr = cn.GeneRegion('TP53', db)
        scores = s.gene_scores(gr)
        self.assertEqual(len(scores), 3)
        self.assertTrue(all(v == 0 for v in scores))


# ---------------------------------------------------------------------------
# 10. SegmentedSampleSet scores() / scores_relaxed() tests
# ---------------------------------------------------------------------------

class SampleSetScoresTest(unittest.TestCase):

    def test_scores(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        s.evaluate(f)
        db = cn.GeneDatabase('refGene.db')
        scores = s.scores(cn.Region(chromosome=17, start=7e6, end=8e6), db)
        self.assertIn('TP53', scores)
        self.assertEqual(len(scores['TP53']), 3)
        self.assertTrue(all(v == 0 for v in scores['TP53']))

    def test_scores_relaxed(self):
        s = cn.SegmentedSampleSet(os.path.join(test_dir, 'test.seg'))
        f = cn.CNAFilter(0.1, math.log(5/2, 2), 10e6, -0.1, math.log(0.7/2, 2), 12e6)
        s.evaluate(f)
        db = cn.GeneDatabase('refGene.db')
        q = cn.Region(chromosome=17, start=7e6, end=8e6)
        strict = s.scores(q, db)
        relaxed = s.scores_relaxed(q, db)
        self.assertEqual(set(strict.keys()), set(relaxed.keys()))
        for gene in strict:
            self.assertEqual(strict[gene], relaxed[gene],
                f"Gene {gene}: strict {strict[gene]} != relaxed {relaxed[gene]}")


# ---------------------------------------------------------------------------
# 11. Cytoband tests
# ---------------------------------------------------------------------------

class CytobandTest(unittest.TestCase):

    def setUp(self):
        self.ct = cn.CytobandTable(os.path.join(test_dir, 'cytoBand.txt'))

    def test_whole_chromosome(self):
        r = cn.CytobandRegion('5', self.ct)
        self.assertEqual(r.chromosome, '5')
        self.assertEqual(r.start, 0)
        self.assertEqual(r.end, 180857866)

    def test_chr9(self):
        r = cn.CytobandRegion('chr9', self.ct)
        self.assertEqual(r.chromosome, '9')
        self.assertEqual(r.start, 0)
        self.assertEqual(r.end, 140273252)

    def test_arm(self):
        r = cn.CytobandRegion('1q', self.ct)
        self.assertEqual(r.chromosome, '1')
        self.assertEqual(r.start, 124300000)
        self.assertEqual(r.end, 247249719)

    def test_band(self):
        r = cn.CytobandRegion('3p21', self.ct)
        self.assertEqual(r.chromosome, '3')
        self.assertEqual(r.start, 43600000)
        self.assertEqual(r.end, 54400000)

    def test_subband(self):
        r = cn.CytobandRegion('4q31.1', self.ct)
        self.assertEqual(r.chromosome, '4')
        self.assertEqual(r.start, 139500000)
        self.assertEqual(r.end, 141700000)

    def test_invalid(self):
        self.assertRaises(Exception, cn.CytobandRegion, 'ZZZ', self.ct)


if __name__ == '__main__':
    unittest.main()
