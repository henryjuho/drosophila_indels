#!/usr/bin/env python

import argparse
from qsub import *

parser = argparse.ArgumentParser()
parser.add_argument('-vcf', '--vcf', help='VCF with variants to filter', required=True)
parser.add_argument('-bed', '--bed_repeats', help='BED file with repeat regions listed', required=True)
parser.add_argument('-ref', '--reference', help='Reference genome', required=True)
parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()

# files
variants = args.vcf
repeats = args.bed_repeats
reference = args.reference
output_prefix = variants.rstrip('.vcf')

# VariantFiltration
VarFil_cmdline = ('java -Xmx6g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
                  '-T VariantFiltration '
                  '-R ' + reference + ' '
                  '-V ' + variants + ' '
                  '-o ' + output_prefix + '.rmarked.vcf '
                  '--mask ' + repeats + ' '
                  '--maskName Repeats')

# SelectVariants
SelVar_cmdline = ('java -Xmx6g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
                  '-T SelectVariants '
                  '-R ' + reference + ' '
                  '-V ' + output_prefix + '.rmarked.vcf '
                  '-o ' + output_prefix + '.rfiltered.pass.vcf '
                  '--excludeFiltered')

# remove intermediate file
rm_cmdline = 'rm ' + output_prefix + '.rmarked.vcf'

# submit job
q_sub([VarFil_cmdline, SelVar_cmdline, rm_cmdline],
      out=output_prefix + '.repeatfiltering',
      mem=15, rmem=10, evolgen=args.evolgen)
