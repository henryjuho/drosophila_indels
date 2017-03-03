#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-bam_list', help='List of bam files to genotype', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
parser.add_argument('-out', help='Output_directory and prefix', required=True)
parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()

# variables
bams = args.bam_list
ref_genome = args.ref
output = args.out + '.samtools.raw.snps.indels.vcf'
evolgen = args.evolgen

# SAMtools
sam_commandline = ('samtools mpileup '
                   '-b ' + bams + ' '
                   '-C 50 '
                   '-f ' + ref_genome + ' '
                   '-u '
                   '| bcftools_new call '
                   '-O v '
                   '-m '
                   '-o ' + output)

q_sub([sam_commandline], out=output, t=168, evolgen=evolgen)
