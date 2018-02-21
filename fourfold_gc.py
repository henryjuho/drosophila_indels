#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
import gzip
import subprocess


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-fourfold', help='bed file of fourfold sites with final column listing gene name',
                        required=True)
    parser.add_argument('-ref', help='Reference genome', required=True)
    parser.add_argument('-gc_thres', help=argparse.SUPPRESS, default=72)
    args = parser.parse_args()

    ref_genome = pysam.FastaFile(args.ref)
    trans_data = pysam.TabixFile(args.fourfold)
    chromosomes = trans_data.contigs

    out_stem = args.fourfold.replace('.bed.gz', '')
    bed_out = '{}_maxgc{}.bed'.format(out_stem, args.gc_thres)
    gc_out = '{}_gc_content.txt'.format(out_stem)

    gene_dict = {}

    # loop through all chromosomes in the fourfold bed
    for chromo in chromosomes:

        ref_str = ref_genome.fetch(chromo)

        # process each bed chromosome
        for line in trans_data.fetch(chromo, parser=pysam.asTuple()):

            chromo, start, stop, trans_id = line[0], int(line[1]), int(line[2]), line[3]

            # add relevant keys
            if chromo not in gene_dict.keys():
                gene_dict[chromo] = {trans_id: [0, 0, 0]}
            if trans_id not in gene_dict[chromo].keys():
                gene_dict[chromo][trans_id] = [0, 0, 0]

            # get ref string for region
            ref_seq = ref_str[start: stop].upper()

            at = ref_seq.count('A') + ref_seq.count('T')
            gc = ref_seq.count('G') + ref_seq.count('C')

            gene_dict[chromo][trans_id][0] += gc
            gene_dict[chromo][trans_id][1] += at
            percent_gc = (gene_dict[chromo][trans_id][0] /
                          float(gene_dict[chromo][trans_id][0] + gene_dict[chromo][trans_id][1])) * 100.0
            gene_dict[chromo][trans_id][2] = percent_gc

    trans_data.close()

    # process gc content
    with open(gc_out, 'w') as gc_file:
        failing_trans = []
        print('transcript\tgc\tat\tpercent_gc')
        for x in gene_dict.keys():
            for transcript in gene_dict[x].keys():
                at_cont, gc_cont, gc_percent = gene_dict[x][transcript]
                print(transcript, at_cont, gc_cont, gc_percent, sep='\t', file=gc_file)

                if gc_percent > args.gc_thres:
                    failing_trans.append(transcript)

    # filter bed
    with open(bed_out, 'w') as bed_file:
        for bed_line in gzip.open(args.fourfold):
            trans_id = bed_line.rstrip().split()[-1]
            if trans_id in failing_trans:
                continue
            else:
                print(bed_line.rstrip(), file=bed_file)

    # bgzip and tabix
    subprocess.call('bgzip {}'.format(bed_out), shell=True)
    subprocess.call('tabix -pbed {}.gz'.format(bed_out), shell=True)


if __name__ == '__main__':
    main()
