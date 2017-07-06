#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import *
import sys
import gzip
import pysam


def get_complement(dna_seq):

    dna_dict = {'A': 'T', 'a': 't', 'T': 'A', 't': 'a',
                'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c',
                'N': 'N', 'n': 'n', '-': '-', '?': '?'}

    comp = ''.join([dna_dict[x] for x in dna_seq])

    return comp


def cds_coord_to_codon_coords(fa_coord, direction):

    iter_coord = [x.split('..') for x in fa_coord.split(',')]
    all_positions = []
    for block in iter_coord:
        block_poss = range(int(block[0]), int(block[1])+1)
        all_positions += block_poss
    if direction == 'complement':
        all_positions = all_positions[::-1]
    return all_positions


def start_stop_ok(fa_seq):

    codon_list = [fa_seq[i:i+3] for i in range(0, len(fa_seq), 3)]
    stops = {'TAG', 'TGA', 'TAA'}
    # check if has start codon
    if codon_list[0] != 'ATG':
        return False

    # check if seq ends in stop codon
    elif codon_list[-1] not in stops:
        return False

    # check for premature stop
    else:
        for c in codon_list[1:-1]:
            if c in stops:
                return False

        return True


def main():
    # arguments
    parser = argparse.ArgumentParser(description='Script that loops through a fasta file of CDS sequence and extracts '
                                                 'alignments for each sequence from a whole genome alignment bed file '
                                                 'and writes each alignment to a separate phylip file.')
    parser.add_argument('-wga_bed',
                        help='Whole genome alignment bed file to extract CDS alignments from',
                        required=True)
    parser.add_argument('-cds_fa',
                        help='Fasta file of CDS transcript sequences',
                        required=True)
    parser.add_argument('-out_dir',
                        help='Output directory to create and write phylip files to',
                        required=True)
    parser.add_argument('-sub',
                        help='If specified will submit script to cluster',
                        action='store_true',
                        default=False)
    parser.add_argument('-evolgen',
                        help='If specified will submit script to lab queue',
                        default=False,
                        action='store_true')
    args = parser.parse_args()

    # submission loop
    if args.sub is True:
        command_line = [' '.join([x for x in sys.argv if x != '-sub' and x != '-evolgen'])]
        q_sub(command_line, out=args.out_dir + 'cds_alignments_from_wga', evolgen=args.evolgen, t=48)
        sys.exit()

    # variables
    fa = args.cds_fa
    wga = pysam.TabixFile(args.wga_bed)
    out_dir = args.out_dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # loop through fasta
    sequence, chromo, coords, skip, method, trans_name = '', '', [], False, '', ''
    for line in gzip.open(fa):

        # skip sequence line following skipped header
        if skip is True:
            skip = False
            continue

        # if a header line
        if line.startswith('>'):

            # skips records with both and forward and reverse strand in coords
            if 'complement' in line and 'join' in line:
                skip = True
                continue

            # process prev sequence
            if sequence != '' and len(sequence) % 3 == 0 and start_stop_ok(sequence):

                align_dict = {}

                # extract positions from alignment
                for codon in [coords[i: i + 3] for i in range(0, len(coords), 3)]:
                    codon_dict = {}
                    for pos in codon:
                        align_pos = list(wga.fetch(chromo, pos - 1, pos, parser=pysam.asTuple()))

                        # get align data in form [('dmel', 'T'), ('dsim', '?'), ('dyak', 'T')]
                        spp_bases = zip(align_pos[0][4].split(','), align_pos[0][7].split(','))
                        for spp in spp_bases:
                            if spp[0] not in codon_dict.keys():
                                codon_dict[spp[0]] = ''
                            if len(align_pos) == 0:
                                codon_dict[spp[0]] += 'N'
                            else:
                                codon_dict[spp[0]] += spp[1][0]

                    # check codons and add to sequence string
                    codon_base_content = ''.join(codon_dict.values())
                    if '-' not in codon_base_content and 'N' not in codon_base_content:
                        for s in codon_dict.keys():
                            if s not in align_dict.keys():
                                align_dict[s] = ''
                            align_dict[s] += codon_dict[s]

                # complement (will have already been reversed) if necessary
                if method == 'complement':
                    align_dict = {x: get_complement(align_dict[x]) for x in align_dict.keys()}

                # output seq file
                seq_len = len(align_dict.values()[0])
                n_spp = len(align_dict.keys())
                file_name = '{}{}.{}.{}'.format(out_dir, trans_name, seq_len, 'phy')
                with open(file_name, 'w') as phy:
                    print('\t{}\t{}'.format(n_spp, seq_len), file_name=phy)
                    for spp_seq in align_dict.keys():
                        out_seq = align_dict[spp_seq]
                        line_wrap = 60
                        print(spp_seq, file_name=phy)
                        for i in range(0, len(out_seq), line_wrap):
                            print(out_seq[i: i + line_wrap], file_name=phy)

            # reset holders and store details of next sequence
            sequence = ''
            header_info = line.split(';')[1]
            chromo = header_info.split(':')[0].replace('loc=', '').replace(' ', '')
            trans_name = line.split(' ')[0].strip('>')
            if '(' in header_info:
                method = header_info.split('(')[0].split(':')[1]
                coords = cds_coord_to_codon_coords(header_info.split('(')[1].rstrip(')'), method)
            else:
                method = 'join'
                coords = cds_coord_to_codon_coords(header_info.split(':')[1], method)

        # if a sequence line add to seq string
        else:
            sequence += line.rstrip()

    # process last seq in file


if __name__ == '__main__':
    main()
