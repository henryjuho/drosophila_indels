#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='vcf file to get mean depth for', required=True)
args = parser.parse_args()

# variables
vcf = args.vcf

# vcftools cmd
vcftools_cmd = 'vcftools --vcf ' + vcf + ' --out depth_stats --depth > ' + vcf.replace('.vcf', '.depth.txt')
q_sub([vcftools_cmd], out=vcf.replace('.vcf', '.depth'))
