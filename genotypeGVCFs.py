#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-ref', help='Reference genome', required=True)
parser.add_argument('-vcf_list', help='List file g.vcfs', required=True)
parser.add_argument('-out', help='Output vcf', required=True)
parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()

# variables
ref = args.ref
vcf_list = args.vcf
out = args.out
evolgen = args.evolgen

# submit job
genotyper = ('java -Xmx6g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
             '--disable_auto_index_creation_and_locking_when_reading_rods '
             '-T GenotypeGVCFs -allSites '
             '-R ' + ref + ' '
             '-V ' + vcf_list + ' '
             '-o ' + out)

q_sub([genotyper], out=out.replace('.vcf', ''), t=60, mem=10, rmem=10, evolgen=evolgen)
