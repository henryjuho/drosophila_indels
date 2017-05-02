#!/usr/bin/env python

from __future__ import print_function
import argparse
import gzip
import sys


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-gff', help='GFF file to read annotations from, note squence names must match VCF',
                        required=True)
    parser.add_argument('-vcf', help='VCF file to annotate variants in', required=True)
    parser.add_argument('-chr', help='Chromosome to annotate, must be consistent with GFF and VCF', required=True)
    args = parser.parse_args()

    # variables
    gff = args.gff
    chromo = args.chr
    vcf_file = args.vcf
    out_vcf = vcf_file.replace('.vcf', '.annotated.' + chromo + '.vcf')

    # make coord sets
    introns = []
    cds = []
    genes = []
    for line in gzip.open(gff):
        line = line.split('\t')
        if line[0] == chromo:
            gff_chromo, feature, feat_start, feat_end = line[0], line[2], int(line[3]), int(line[4])

            if chromo == gff_chromo:
                if feature == 'gene':
                    genes += range(feat_start-1, feat_end)
                elif feature == 'intron':
                    introns += range(feat_start-1, feat_end)
                elif feature == 'CDS':
                    cds += range(feat_start-1, feat_end)
                else:
                    continue

    introns = set(introns)
    cds = set(cds)
    genes = set(genes)

    # annotation counters
    all_variants = 0
    cds_frame = 0
    cds_nonframe = 0
    inter = 0
    intron_count = 0

    # loop through vcf and identify category for each variant
    previous_line = ''
    with open(out_vcf, 'w') as annotated_vcf:
        for line in open(vcf_file):
            if line.startswith('#'):
                if line.startswith('##contig') and previous_line.startswith('##INFO'):
                    new_info = '##INFO=<ID=ANNO,Number=1,Type=String,Description="Annotation of genomic region">\n'
                    annotated_vcf.write(new_info)
                previous_line = line
                annotated_vcf.write(line)
            elif line.split('\t')[0] == chromo:
                all_variants += 1
                split_line = line.split('\t')
                chromo, start, end = split_line[0], int(split_line[1])-1, int(split_line[1]) + len(split_line[3])

                in_cds = False
                for x in range(start, end):
                    if x not in cds:
                        break
                    else:
                        in_cds = True

                in_intergenic = False
                for x in range(start, end):
                    if x in genes:
                        break
                    else:
                        in_intergenic = True

                in_intron = False
                for x in range(start, end):
                    if x not in introns:
                        break
                    else:
                        in_intron = True

                region_true = [x for x in [in_cds, in_intergenic, in_intron] if x is True]
                number_true = len(region_true)

                if number_true == 1:
                    if in_cds is True:
                        variant_length = abs(len(split_line[3]) - len(split_line[4]))
                        if variant_length % 3.0 == 0.0:
                            region = 'CDS_non_frameshift'
                            cds_nonframe += 1
                        else:
                            region = 'CDS_frameshift'
                            cds_frame += 1
                    elif in_intron is True:
                        intron_count += 1
                        region = 'intron'
                    elif in_intergenic is True:
                        inter += 1
                        region = 'intergenic'
                    else:
                        print('BUGGY BUGGY BUGGY!')
                        sys.exit('An impossible error occured!')

                    # write newly annotated line
                    anno_line = ('\t'.join(split_line[0:7]) + '\t' + split_line[7] +
                                 ';ANNO=' + region + '\t' + '\t'.join(split_line[8:]))
                    annotated_vcf.write(anno_line)
                else:
                    annotated_vcf.write(line)
            else:
                continue

    print('\n'
          '|Category          | Number INDELs   |\n'
          '|:-----------------|:---------------:|\n'
          '|All               | ' + str(all_variants) + ' |\n'
          '|CDS_frameshift    | ' + str(cds_frame) + ' |\n'
          '|CDS_non_frameshift| ' + str(cds_nonframe) + ' |\n'
          '|Intron            | ' + str(intron_count) + ' |\n'
          '|Intergenic        | ' + str(inter) + ' |\n'
          '|Not annotated     | ' + str(all_variants - (cds_frame + cds_nonframe + intron_count + inter)) + ' |')

if __name__ == '__main__':
    main()
