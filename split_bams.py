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

# split each bam
counter = 0
for bam in bams:

    bam_basename = bam[bam.rfind('/')+1:].replace('.bam', '')
    get_id_cmd = 'samtools view -H ' + bam + ' | grep @RG | cut -f 6 | uniq | cut -d ":" -f 2'
    ids = subprocess.Popen(get_id_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]

    for indiv in ids:
        rg_file = out_dir + indiv + '.' + bam_basename + '.rg.txt'
        new_bam = out_dir + indiv + '.' + bam_basename + '.bam'
        rg_id_grep = 'samtools view -H ' + bam + ' | grep @RG | grep SM:' + indiv + ' | cut -f 2 | cut -d ":" -f 2'
        rgs = subprocess.Popen(rg_id_grep, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
        with open(rg_file, 'w') as rg_txt:
            print('\n'.join(rgs), file=rg_txt)

        split_cmd = 'samtools view -bhR ' + rg_file + ' ' + bam + ' > ' + new_bam

        q_sub([split_cmd], out=out_dir + indiv + '.' + bam_basename, evolgen=evolgen)
