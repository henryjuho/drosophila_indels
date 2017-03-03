#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-ref', help='Reference genome', required=True)
parser.add_argument('-bam_dir', help='Directory containing BAM files', required=True)
parser.add_argument('-out_dir', help='Output directory', required=True)
parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()

# variables
ref = args.ref
bam_dir = args.bam_dir
bam_data = [(bam_dir + x, x.split('.')[0]) for x in os.listdir(bam_dir) if x.endswith('.bam')]
out_dir = args.out_dir
if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
evolgen = args.evolgen

# per bam command
for info in bam_data:
    bam = info[0]
    sample = info[1]
    hap_caller = ('java -Xmx4g -jar ~/gatk3.7/GenomeAnalysisTK.jar '
                  '-T HaplotypeCaller '
                  '-R ' + ref + ' '
                  '-I ' + bam + ' '
                  '--emitRefConfidence GVCF '
                  '-ploidy 2 '
                  '-o ' + out_dir + sample + '.raw.snps.indels.g.vcf')

    q_print([hap_caller], out=out_dir + sample + '.hap_calling', t=24, evolgen=evolgen)
