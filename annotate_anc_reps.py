#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import *
import sys
import gzip


def in_region(reg_start, reg_end, set_name):
    for x in range(reg_start, reg_end):
        if x not in set_name:
            return False
        else:
            continue
    return True


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', help='bed file containing ancestral repeats coordinates', required=True)
    parser.add_argument('-vcf', help='vcf file to annotate', required=True)
    parser.add_argument('-trim_non_anc_reps',
                        help='If specified will remove variants that are in non ancestral repeats',
                        action='store_true', default=False)
    parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
    args = parser.parse_args()

    # submission loop
    if args.sub is True:
        command_line = [' '.join([x for x in sys.argv if x != '-sub'])]
        q_sub(command_line, out=args.vcf+'ar')
        sys.exit(0)

    # variables
    bed = args.bed
    vcf = args.vcf
    anno_vcf = vcf.replace('.vcf', '.ar.vcf')
    trim = args.trim_non_anc_reps

    # repeats dictionary
    rep_dict = {}
    for line in gzip.open(bed):
        chromo, start, end = line.split('\t')[0], int(line.split('\t')[1]), int(line.split('\t')[2])
        if chromo not in rep_dict.keys():
            rep_dict[chromo] = range(start, end)
        else:
            rep_dict[chromo] += range(start, end)
    rep_dict = {x[0]: set(x[1]) for x in rep_dict.items()}

    # adjust anno in vcf
    with open(anno_vcf, 'w') as out_vcf:
        for line in open(vcf):
            if line.startswith('#'):
                print(line, file=out_vcf)
            else:
                vcf_chromo, vcf_start, vcf_var_len = line.split('\t')[0], int(line.split('\t')[1]), \
                                                     len(line.split('\t')[3])
                vcf_end = vcf_start + vcf_var_len
                if in_region(vcf_start-1, vcf_end, rep_dict[vcf_chromo]) is False:
                    if trim is True and 'Repeats' in line:
                        continue
                    print(line, file=out_vcf)
                else:
                    if 'ANNO=intergenic' not in line:
                        if trim is True:
                            continue
                        print(line, file=out_vcf)
                    else:
                        line = line.replace('ANNO=intergenic', 'ANNO=intergenic_ar')
                        print(line, file=out_vcf)

if __name__ == '__main__':
    main()
