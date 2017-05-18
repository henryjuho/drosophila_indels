#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
import gzip
import itertools


def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
        b = list(b)
        yield b[0][1], b[-1][1]


def gff_regions(gff_file, chromo):

    # dict holding region coords for each chromosome
    features = ['CDS', 'intron', 'gene', 'all_feat']
    region_dict = {y: [] for y in features}

    # make coord sets
    for line in gzip.open(gff_file):
        line = line.split('\t')
        if line[0] == chromo:
            gff_chromo, feature, feat_start, feat_end = line[0], line[2], int(line[3]), int(line[4])

            region_dict['all_feat'] += range(feat_start - 1, feat_end)
            if feature in features:
                region_dict[feature] += range(feat_start - 1, feat_end)

    region_dict = {k2: set(region_dict[k2]) for k2 in region_dict.keys()}

    return region_dict


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-gff', help='gff file to get genomic region coords from', required=True)
    parser.add_argument('-call_fa', help='Callable fasta file', required=True)
    parser.add_argument('-chr_list', help='file of chromosomes to calc callable sites for', required=True)
    args = parser.parse_args()

    # variables
    gff = args. gff
    fa = pysam.FastaFile(args.call_fa)
    chr_list = [x.rstrip() for x in open(args.chr_list)]
    regions = ['all', 'CDS', 'intron', 'intergenic', 'AR']
    call_data = {x: {y: 0 for y in regions} for x in ['total'] + chr_list}  # {chromo: {all: 0, CDS: 0, intron: 0 ...}}

    # get call sites for all chr and regions
    for chromo in chr_list:
        fasta_string = fa.fetch(chromo)
        region_coords = gff_regions(gff, chromo)
        for region in regions:
            if region == 'all':
                callable_sites = fasta_string.count('2')
            elif region == 'AR':
                callable_sites = fasta_string.count('3')
            elif region == 'intergenic':

                # need to reverse engineer intergenic coords
                gene_coords = region_coords['gene']
                if len(gene_coords) == 0:
                    gene_coords = region_coords['all_feat']

                # invert
                intergenic_coords = [i for i in range(0, len(fasta_string)) if i not in gene_coords]
                intergenic_coords = list(ranges(intergenic_coords))

                callable_sites = ''.join([fasta_string[x[0]: x[1]+1] for x in intergenic_coords]).count('2')

            else:
                coords = list(ranges(sorted(list(region_coords[chromo][region]))))
                callable_sites = ''.join([fasta_string[x[0]: x[1] + 1] for x in coords]).count('2')

            call_data[chromo][region] += callable_sites
            call_data['total'][region] += callable_sites

    # output sites
    print('\t'.join(['contig', 'region', 'no_callable']))
    for seq in sorted(call_data.keys()):
        for reg in call_data[seq].keys():
            print('\t'.join([seq, reg, str(call_data[seq][reg])]))

if __name__ == '__main__':
    main()
