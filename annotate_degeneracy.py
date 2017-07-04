#!/usr/bin/env python

from __future__ import print_function
import argparse
import gzip


def bed_to_region_dict(bed):

    coord_dict = {}

    for line in gzip.open(bed):
        line = line.split()
        chromo, start, stop = line[0], int(line[1]), int(line[2])
        if chromo not in coord_dict.keys():
            coord_dict[chromo] = set()
        coord_dict[chromo] |= set(range(start+1, stop+1))

    return coord_dict


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='VCF file to annotate variants in', required=True)
    parser.add_argument('-zerofold', help='bed file of zerofold sites', required=True)
    parser.add_argument('-fourfold', help='bed file of fourfold sites', required=True)
    args = parser.parse_args()

    # variables
    vcf = args.vcf
    zerofold_sites = bed_to_region_dict(args.zerofold)
    fourfold_sites = bed_to_region_dict(args.fourfold)
    annotated_vcf = vcf.replace('.vcf', '.degen.vcf')

    # loop through vcf
    zero_count = 0
    four_count = 0

    with open(annotated_vcf, 'w') as out_vcf:
        previous_line = ''
        for line in open(vcf):
            if line.startswith('#'):
                if line.startswith('##contig') and previous_line.startswith('##INFO'):
                    new_info = '##INFO=<ID=DEGEN,Number=1,Type=Integer,Description="Annotation of SNP degeneracy">'
                    print(new_info, file=out_vcf)
                previous_line = line
                print(line.rstrip(), file=out_vcf)
            else:
                split_line = line.rstrip().split('\t')
                chromo = split_line[0]
                snp_pos = int(split_line[1])

                # sets line unchanged if chromo doesn't appear in bed files
                if chromo not in zerofold_sites.keys():
                    out_line = line.rstrip()

                # if snp is zerofold sets new line
                elif snp_pos in zerofold_sites[chromo]:
                    out_line = '\t'.join(split_line[0:7] + [split_line[7] + ';DEGEN=0'] + split_line[8:])
                    zero_count += 1

                # if snp is fourfold sets new line
                elif snp_pos in fourfold_sites[chromo]:
                    out_line = '\t'.join(split_line[0:7] + [split_line[7] + ';DEGEN=4'] + split_line[8:])
                    four_count += 1

                # if snp is neither sets original line
                else:
                    out_line = line.rstrip()

                # outputs line
                print(out_line, file=out_vcf)

    summary_table = ('\n'
                     '|Category          | Number SNPs     |\n'
                     '|:-----------------|:---------------:|\n'
                     '|zerofold          |{:<17}|\n'
                     '|fourfold          |{:<17}|').format(zero_count, four_count)

    print(summary_table)

if __name__ == '__main__':
    main()
