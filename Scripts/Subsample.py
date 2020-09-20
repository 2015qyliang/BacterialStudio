# coding: utf-8

from Bio import SeqIO
import random
import sys

subsampleSize = sys.argv[1]
inputFile = sys.argv[2]
outputFile = sys.argv[3]

random.seed(123)

fileSeqs = { str(seq.id):str(seq.seq)  for seq in SeqIO.parse(inputFile, 'fasta') }

subsampleID = random.sample(list(fileSeqs.keys()), int(subsampleSize))

newSeqs = [ '>' + str(seqid) + '\n' + str(fileSeqs[seqid]) + '\n' for seqid in subsampleID ]

open(outputFile, 'w').writelines(newSeqs)
