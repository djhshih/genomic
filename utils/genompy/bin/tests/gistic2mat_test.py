#!/usr/bin/env python3

import math, subprocess

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


def testMain():
	subprocess.call(
		['../gistic2mat.py', 'test3.txt', 'test3.cn',
			'--segfile', 'test3.seg', '--genesDb', 'refGene.db', '--genefile', 'test3.genes'] )
	checkFiles('test3.cn', 'test3.cn.ans')
	checkFiles('test3.genes', 'test3.genes.ans')

if __name__ == '__main__':
	testMain()

