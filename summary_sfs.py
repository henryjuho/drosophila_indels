#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='vcf file to summarise sfs for', required=True)
parser.add_argument('-mode', help='mode to run in', choices=['SNP', 'INDEL'], required=True)
parser.add_argument('-out', help='Output directory', required=True)
parser.add_argument('-spp', help='Species', required=True)
args = parser.parse_args()

# variables
vcf = args.vcf
out = args.out
spp = args.spp
mode = args.mode
mode_dict = {'SNP': ['snp'], 'INDEL': ['ins', 'del']}
region_dict = {'SNP': ['ALL', 'CDS_non_frameshift', 'intron', 'intergenic', 'intergenic_ar'],
               'INDEL': ['ALL', 'CDS_frameshift', 'CDS_non_frameshift', 'intron', 'intergenic', 'intergenic_ar']}
# SFS
for m in mode_dict[mode]:
    for r in region_dict[mode]:
        new_out = '{}{}_{}_{}.sfs.txt'.format(out, spp, m, r.replace('_', ''))
        sfs_cmd = ('~/sfs_utils/vcf2raw_sfs.py '
                   '-vcf {} '
                   '-mode {} '
                   '-auto_only '
                   '-region {} '
                   '| sort | uniq -c | while read i; do echo " "$i; done | tr -s " " '
                   '> {}').format(vcf, m, r, new_out)
        q_sub([sfs_cmd], out=new_out.replace('.txt', ''))
