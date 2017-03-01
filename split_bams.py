#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-bam_dir', help='directory containing bam files', required=True)
parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()

# variables
bam_dir = args.bam_dir
bams = [bam_dir + x for x in os.listdir(bam_dir) if x.endswith('.bam')]
evolgen = args.evolgen

# split each bam
counter = 0
for bam in bams:
    counter += 1
    split_cmd = "samtools split -f '%!.%*.bam' " + bam
    q_sub([split_cmd], out=bam_dir + 'bam_split_' + str(counter), evolgen=evolgen)
