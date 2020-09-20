# coding: utf-8

import sys
import os

inputFolder = sys.argv[1]
outputFolder = sys.argv[2]

open(outputFolder + 'RST6_AllNochimeras.fasta', 'w')

for fn in os.listdir(inputFolder):
	seqs = [ line for line in open(inputFolder + fn).readlines() ]
	open(outputFolder + 'RST6_AllNochimeras.fasta', 'a').writelines(seqs)
