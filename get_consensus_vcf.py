#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf_I', help='First VCF to compare, and to use for naming output', required=True)
parser.add_argument('-vcf_II', help='Second VCF to compare', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
parser.add_argument('-out', help='Output_directory', required=True)
parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()

# variables
vcf_1 = args.vcf_I
vcf_2 = args.vcf_II
ref_genome = args.ref
out = args.out
evolgen = args.evolgen

for var in ['SNP', 'INDEL']:

    # GATK select variants to extract variants
    vcf_1_out = args.out + vcf_1[vcf_1.rfind('/')+1:].replace('.allsites.vcf', '.raw.' + var.lower() + 's.vcf')
    extract1_cmd = ('java -Xmx6g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
                    '-T SelectVariants '
                    '-R ' + ref_genome + ' '
                    '-V ' + vcf_1 + ' '
                    '-selectType ' + var + ' '
                    '-trimAlternates '
                    '-env '
                    '-o ' + vcf_1_out)

    vcf_2_out = args.out + vcf_2[vcf_2.rfind('/')+1:].replace('.allsites.vcf', '.raw.' + var.lower() + 's.vcf')
    extract2_cmd = ('java -Xmx6g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
                    '-T SelectVariants '
                    '-R ' + ref_genome + ' '
                    '-V ' + vcf_2 + ' '
                    '-selectType ' + var + ' '
                    '-trimAlternates '
                    '-env '
                    '-o ' + vcf_2_out)

    # GATK select variants with concordance set
    consensus_file = args.out + vcf_1[vcf_1.rfind('/')+1:].split('.')[0] + '.consensus.raw.' + var.lower() + 's.vcf'
    concensus_cmd = ('java -Xmx6g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
                     '-T SelectVariants '
                     '-R ' + ref_genome + ' '
                     '-V ' + vcf_1_out + ' '
                     '--concordance ' + vcf_2_out + ' '
                     '-o ' + consensus_file)

    # submit to cluster
    q_sub([extract1_cmd, extract2_cmd, concensus_cmd],
          out=args.out + vcf_1[vcf_1.rfind('/')+1:].split('.')[0] + '.consensus.raw' + var.lower() + 's',
          mem=10, rmem=10, t=72, evolgen=evolgen)
