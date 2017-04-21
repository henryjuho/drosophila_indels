#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='vcf file to exclude SNPs overlapping INDELs', required=True)
args = parser.parse_args()

out_vcf = args.vcf.replace('.vcf', '.exsnpindel.vcf')

# the power of grep
grep_1 = 'grep ^# {} > {}'.format(args.vcf, out_vcf)
grep_2 = 'grep -v ^# {} | grep -v "*" >> {}'.format(args.vcf, out_vcf)

q_sub([grep_1, grep_2], out=args.vcf.replace('vcf', 'exclude_snp_in_indel'))
