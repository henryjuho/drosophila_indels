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
    region_dict = {y: set() for y in features}

    # make coord sets
    for line in gzip.open(gff_file):
        line = line.split('\t')
        if line[0] == chromo:
            gff_chromo, feature, feat_start, feat_end = line[0], line[2], int(line[3]), int(line[4])

            region_dict['all_feat'] |= set(range(feat_start - 1, feat_end))
            if feature in features:
                region_dict[feature] |= set(range(feat_start - 1, feat_end))

    # region_dict = {k2: [x for x in ranges(sorted(list(region_dict[k2])))] for k2 in region_dict.keys()}

    return region_dict


def bed_regions(bed_file, chromo):
    open_bed = pysam.TabixFile(bed_file)
    coords = []
    for line in open_bed.fetch(chromo, parser=pysam.asTuple()):
        start = int(line[1])
        end = int(line[2])
        coords += range(start, end)

    return coords


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-gff', help='Gff file to get genomic region coords from', required=True)
    parser.add_argument('-call_fa', help='Callable fasta file', required=True)
    parser.add_argument('-chr_list', help='File of chromosomes to calc callable sites for', required=True)
    parser.add_argument('-opt_bed', help='Optional bed file of regions to count callables sites, with associated label'
                                         'i.e. /path/to/file.bed.gz,my_sub_region', action='append')
    args = parser.parse_args()

    # variables
    gff = args.gff
    fa = pysam.FastaFile(args.call_fa)
    if args.opt_bed is not None:
        bed_files = {x.split(',')[1]: x.split(',')[0] for x in args.opt_bed}
    else:
        bed_files = {}

    chr_list = [x.rstrip() for x in open(args.chr_list)]
    regions = ['ALL', 'CDS', 'intron', 'intergenic', 'AR'] + bed_files.keys()
    call_data = {x: {y: {'ALL': 0, 'POL': 0} for y in regions} for x in ['ALL'] + chr_list}

    # {chromo: {all: 0, CDS: 0, intron: 0 ...}}

    # get call sites for all chr and regions
    for chromo in chr_list:
        fasta_string = fa.fetch(chromo)
        region_coords = gff_regions(gff, chromo)
        for region in regions:
            if region == 'ALL':
                callable_sites_all = fasta_string.upper().count('K')
                callable_sites_pol = fasta_string.count('K')
            elif region == 'AR':
                callable_sites_all = fasta_string.upper().count('R')
                callable_sites_pol = fasta_string.count('R')
            elif region == 'intergenic':
                # need to reverse engineer intergenic coords
                gene_coords = region_coords['gene']
                if len(gene_coords) == 0:
                    gene_coords = region_coords['all_feat']

                # invert
                intergenic_coords = [i for i in range(0, len(fasta_string)) if i not in gene_coords]
                intergenic_coords = list(intergenic_coords)

                callable_sites_all = ''.join([fasta_string[x] for x in intergenic_coords]).upper().count('K')
                callable_sites_pol = ''.join([fasta_string[x] for x in intergenic_coords]).count('K')

            # handle optional bed files
            elif region in bed_files.keys():
                degen_bed = bed_files[region]
                try:
                    degen_coords = bed_regions(degen_bed, chromo)

                    callable_sites_all = ''.join([fasta_string[x] for x in degen_coords]).upper().count('K')
                    callable_sites_pol = ''.join([fasta_string[x] for x in degen_coords]).count('K')
                except ValueError:
                    callable_sites_all = 0
                    callable_sites_pol = 0

            else:
                coords = list(region_coords[region])
                callable_sites_all = ''.join([fasta_string[x] for x in coords]).upper().count('K')
                callable_sites_pol = ''.join([fasta_string[x] for x in coords]).count('K')

            call_data[chromo][region]['ALL'] += callable_sites_all
            call_data[chromo][region]['POL'] += callable_sites_pol

            call_data['ALL'][region]['ALL'] += callable_sites_all
            call_data['ALL'][region]['POL'] += callable_sites_pol

    # output sites
    print(','.join(['contig', 'region', 'all_callable', 'pol_callable']))
    for seq in sorted(call_data.keys()):
        for reg in call_data[seq].keys():
            print(','.join([seq, reg, str(call_data[seq][reg]['ALL']), str(call_data[seq][reg]['POL'])]))

if __name__ == '__main__':
    main()
