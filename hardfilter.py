#!/usr/bin/env python

import argparse
from qsub import *

# command line options
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', '--VCF_input', help='VCF file to filter', required=True)
parser.add_argument('-ref', '--reference', help='Reference genome', required=True)
parser.add_argument('-mode', help='Mode to run in, SNP or INDEL', choices=['SNP', 'INDEL'], required=True)
parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()


# variables
input_vcf = args.VCF_input
output_prefix = ''
if input_vcf.endswith('.vcf.gz'):
    output_prefix = input_vcf.rstrip('.vcf.gz')
elif input_vcf.endswith('.vcf'):
    output_prefix = input_vcf.rstrip('.vcf')
ref = args.reference
mode = args.mode

# GATK command lines
hard_filter = ('java -Xmx6g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
               '-T VariantFiltration '
               '-R ' + ref + ' '
               '-V ' + input_vcf + ' '
               '-o ' + output_prefix + '.hmarked.vcf ')
if mode == 'INDEL':
    hard_filter += ('--filterExpression "QD<2.0" --filterName "indelQD" '
                    '--filterExpression "FS>200.0" --filterName "indelFS" '
                    '--filterExpression "ReadPosRankSum<-20.0" --filterName "indelRPRS" '
                    '--filterExpression "SOR>10.0" --filterName "indelSOR"')

else:
    hard_filter += ('--filterExpression "QD<2.0" --filterName "snpQD" '
                    '--filterExpression "MQ<40.0" --filterName "snpMQ" '
                    '--filterExpression "FS>60.0" --filterName "snpFS" '
                    '--filterExpression "MQRankSum<-12.5" --filterName "snpMQRS" '
                    '--filterExpression "ReadPosRankSum<-8.0" --filterName "snpRPRS" '
                    '--filterExpression "SOR>3.0" --filterName "snpSOR"')

extract_passed = ('java -Xmx6g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
                  '-T SelectVariants '
                  '-R ' + ref + ' '
                  '-V ' + output_prefix + '.hmarked.vcf '
                  '-o ' + output_prefix + '.hfiltered.vcf '
                  '--excludeFiltered')

# write qsub job
q_sub([hard_filter, extract_passed], out=output_prefix + '.filtering', mem=10, rmem=10, evolgen=args.evolgen)
