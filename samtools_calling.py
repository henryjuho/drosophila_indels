#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-bam_list', help='List of bam files to genotype', required=True)
parser.add_argument('-ref', help='Reference genome location', required=True)
parser.add_argument('-out', help='Output_directory and prefix', required=True)
parser.add_argument('-per_chr', help='If specified will split calling per chromosome',
                    default=False, action='store_true')
parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()

# variables
bams = args.bam_list
ref_genome = args.ref
output = args.out + '.samtools.allsites.vcf'
per_chr = args.per_chr
evolgen = args.evolgen

# SAMtools whole genome
if per_chr is False:
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

# SAMtools per chromosome
else:
    # generate chromosome list
    chromosome_scaffold_list = [[line.split()[0], line.split()[1]] for line in open(ref_genome + '.fai')]
    chromosome_list = [entry[0] for entry in chromosome_scaffold_list if not entry[0].startswith('NODE')] + \
                      ['scaffolds']
    scaffold_list = [entry for entry in chromosome_scaffold_list if entry[0].startswith('NODE')]

    # create scaffold bed file
    scaffold_file_name = args.out + '.scaffolds.bed'
    with open(scaffold_file_name, 'w') as scaffold_file:
        for scaffold in scaffold_list:
            scaffold_file.write(scaffold[0] + '\t0\t' + scaffold[1] + '\n')

    # generate samtools job for each chromosome
    job_list = []
    vcf_list = []
    for position in chromosome_list:
        new_output = args.out + '.samtools.allsites.' + position + '.vcf'
        vcf_list.append(new_output)
        job_name = 'samtools.' + position + '.sh'
        job_list.append(job_name)
        if position != 'scaffolds':
            # SAMtools
            SAM_commandline = ('samtools mpileup '
                               '-b ' + bams + ' '
                               '-C 50 '
                               '-f ' + ref_genome + ' '
                               '-r ' + position + ' '
                               '-u '
                               '| bcftools_new call '
                               '-O v '
                               '-m '
                               '-o ' + new_output)

            # submit
            q_sub([SAM_commandline], out=new_output.replace('.vcf', ''), t=168, jid=job_name, evolgen=evolgen)

        else:
            # SAMtools
            SAM_commandline = ('samtools mpileup '
                               '-b ' + bams + ' '
                               '-C 50 '
                               '-f ' + ref_genome + ' '
                               '-l ' + scaffold_file_name + ' '
                               '-u '
                               '| bcftools_new call '
                               '-O v '
                               '-m '
                               '-o ' + new_output)
            # submit
            q_sub([SAM_commandline], out=new_output.replace('.vcf', ''), t=168, jid=job_name, evolgen=evolgen)

    # concat job
    vcf_list_file = args.out + '.vcf_list.txt'
    with open(vcf_list_file, 'w') as vcf_info:
        for vcf in vcf_list:
            vcf_info.write(vcf + '\n')

    concat_cmd = ('bcftools_new concat -f ' + vcf_list_file + ' -O v -o ' + output)
    q_sub([concat_cmd], out=args.out + '.vcf_merging', hold=job_list, evolgen=evolgen)
