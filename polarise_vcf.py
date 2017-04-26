#!/usr/bin/env python

import argparse
from qsub import *
import sys
import pysam


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='VCF file containing insertion and deletion variants to polarise', required=True)
    parser.add_argument('-wga_bed', help='tabix indexed whole genome alignment bed file', required=True)
    parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
    args = parser.parse_args()

    # submission loop
    if args.sub is True:
        command_line = [' '.join([x for x in sys.argv if x != '-sub'])]
        q_sub(command_line, out=args.vcf.replace('.vcf', 'polarisation'))
        sys.exit()

    # variables
    vcf_file = args.vcf
    out_vcf = vcf_file.replace('.vcf', '.polarised.vcf')
    wga_bed = pysam.TabixFile(args.wga_bed)

    # counters
    counter = 0
    match = 0
    no_hotspot = 0
    low_coverage = 0
    ambiguous = 0

    # loop through vcf file and annotate INDELs
    previous_line = ''
    with open(out_vcf, 'w') as annotated_vcf:
        for line in open(vcf_file):

            # added new header line
            if line.startswith('#'):
                if line.startswith('##contig') and previous_line.startswith('##INFO'):
                    new_info = '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\n'
                    annotated_vcf.write(new_info)
                previous_line = line
                annotated_vcf.write(line)

            # process variants
            else:
                counter += 1
                orig_line = line
                line = line.split('\t')
                chrom, pos, ref, alt, info = line[0], line[1], line[3], line[4], line[7]

                if len(ref) > 1:
                    var_type = 'INDEL'
                else:
                    var_type = 'SNP'

                # get aligned seq for var
                var_align = [x for x in wga_bed.fetch(chrom, int(pos) - 1, int(pos) - 1 + len(ref),
                                                      parser=pysam.asTuple())]

                # skip if not in alignment
                if len(var_align) == 0:
                    annotated_vcf.write(orig_line)
                    low_coverage += 1
                    continue

                elif len(var_align) > 1:
                    if var_type == 'SNP':
                        annotated_vcf.write(orig_line)
                        no_hotspot += 1
                        continue

                    # catch deletions rel to ref that uniq to ref spp
                    else:
                        # merge sequences from multiple bed rows
                        merged_align = [''.join(y) for y in zip(*[x[7].split(',') for x in var_align])]
                        if '-' not in ''.join(merged_align):
                            var_align = merged_align
                        else:
                            annotated_vcf.write(orig_line)
                            no_hotspot += 1
                            continue

                else:
                    var_align = [x[7] for x in var_align][0].split(',')

                var_align = [x.upper() for x in var_align]
                # print chrom, pos, len(ref), ref, alt, var_align

                # skip positions without full coverage
                if '?' in ''.join(var_align):
                    annotated_vcf.write(orig_line)
                    low_coverage += 1
                    continue

                # skip indel hotspots
                indel_sequences = [y.rstrip('-') for y in var_align]
                if '-' in ''.join(indel_sequences):
                    annotated_vcf.write(orig_line)
                    no_hotspot += 1
                    continue

                # skips sites where ref allele differs from that in alignment, ie insertion within INDEL
                if ref != var_align[0].upper():
                    no_hotspot += 1
                    annotated_vcf.write(orig_line)
                    continue

                else:
                    # identify if ref or alt is ancestral
                    out_group_seqs = var_align[1:]
                    ref_anc = True
                    alt_anc = True

                    # identify ref ancestral
                    for sequence in out_group_seqs:
                        if var_type == 'INDEL':
                            if len(ref) != len(sequence.rstrip('-')):
                                ref_anc = False
                                break
                        else:
                            if ref != sequence.rstrip('-'):
                                ref_anc = False
                                break

                    # identify alt ancestral
                    for sequence in out_group_seqs:
                        if var_type == 'INDEL':
                            if len(alt) != len(sequence.rstrip('-')):
                                alt_anc = False
                                break
                        else:
                            if alt != sequence.rstrip('-'):
                                alt_anc = False
                                break

                    # skip ambiguous sites
                    if alt_anc is ref_anc:
                        ambiguous += 1
                        annotated_vcf.write(orig_line)
                        continue

                    # set AA annotation
                    aa = 'NONE'
                    if ref_anc is True:
                        aa = ';AA=' + ref
                    elif alt_anc is True:
                        aa = ';AA=' + alt

                    # write annotation
                    if not aa == 'NONE':
                        polarised_line = '\t'.join(line[0:7]) + '\t' + info + aa + '\t' + '\t'.join(line[8:])
                        annotated_vcf.write(polarised_line)
                        match += 1
                        # print ref, alt, aa

    print('\n'
          'Total no INDELs  : ' + str(counter) + '\n'
          'INDELs polarised : ' + str(match) + '\n'
          'Hotspots         : ' + str(no_hotspot) + '\n'
          'Low spp coverage : ' + str(low_coverage) + '\n'
          'Ambiguous        : ' + str(ambiguous) + '\n'
          'Total unpolarised: ' + str(counter - match) + '\n')

if __name__ == '__main__':
    main()
