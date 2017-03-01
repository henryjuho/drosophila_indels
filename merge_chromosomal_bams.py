#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-bam_dir', help='directory containing bam files', required=True)
parser.add_argument('-out_dir', help='output_directory', required=True)
parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()

# variables
bam_dir = args.bam_dir
bams = [bam_dir + x for x in os.listdir(bam_dir) if x.endswith('.bam')]
out_dir = args.out_dir
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
evolgen = args.evolgen

# make id list
ids = set([x.split('.')[0] for x in os.listdir(bam_dir) if x.endswith('.bam')])

# loop through each individual
for fly in ids:
    whole_genome_bam = out_dir + fly + '.realigned.bam'
    list_file = out_dir + fly + '.bam_list.txt'
    chr_list_make = 'ls ' + bam_dir + fly + '*.bam > ' + list_file
    subprocess.call(chr_list_make, shell=True)

    merge_cmd = 'samtools merge -cpb ' + list_file + ' ' + whole_genome_bam
    q_sub([merge_cmd], out=out_dir + fly + '.chrmerge', evolgen=evolgen)
