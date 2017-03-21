#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='VCF to filter', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
args = parser.parse_args()

# variables
input_vcf = args.vcf
output_vcf_prefix = input_vcf.rstrip('vcf')+'50bp_max.bial'
ref_genome = args.ref

# gatk select variants
select_var_cmdline = ('java -Xmx6g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
                      '-T SelectVariants '
                      '-R ' + ref_genome + ' '
                      '-V ' + input_vcf + ' '
                      '-o ' + output_vcf_prefix + '.vcf '
                      '-restrictAllelesTo BIALLELIC '
                      '--maxIndelSize 50')

# submit job
q_sub([select_var_cmdline], out=output_vcf_prefix, rmem=10, mem=10)
