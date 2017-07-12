#!/usr/bin/env python

from __future__ import print_function
import argparse
import os


def read_phylip(phy_file):
    sequence, spp, out_data = '', '', []
    counter = 0

    for line in open(phy_file):
        if line.startswith('\t'):
            n_seq, len_seq = line.split('\t')[1], int(line.split('\t')[2])
            out_data += [n_seq, len_seq]
        elif line[0].upper() in {'A', 'T', 'G', 'C', 'N', '?', '-'}:
            sequence += line.rstrip()
        else:
            counter += 1

            if counter != 1:
                out_data += [(spp, sequence)]
                sequence = ''

            spp = line.rstrip()

    out_data += [(spp, sequence)]

    return out_data


def main():
    # arguments
    parser = argparse.ArgumentParser(description='Script that trims the final codon (assumed to be the stop codon, '
                                                 'and also trims any codons containing "?" or "-". Outputs new phylip '
                                                 'files to specified directory. Note any input files with entire '
                                                 'sequence missing for a species will not be output.')
    parser.add_argument('-phy_dir', help='Directory containg phylip files to trim', required=True)
    parser.add_argument('-out_dir', help='Output directory to write new phylip files to', required=True)
    args = parser.parse_args()

    # loop through alignments
    for phy in [x for x in os.listdir(args.phy_dir) if x.endswith('.phy')]:
        in_phy = args.phy_dir + phy
        out_phy = args.out_dir + phy.replace('.phy', '.trim.phy')

        align_data = read_phylip(in_phy)
        # get seqs into form [(spp1, spp2, spp3), (se1, seq2, seq3)]
        sequences = zip(*align_data[2:])

        trimmed_sequences = ['' for i in range(0, len(sequences[0]))]
        removed_bases = 0

        for i in range(0, len(sequences[1][0]), 3):
            codons = [c[i: i+3] for c in sequences[1]]

            # check codons - remove if contain ? or -
            if '-' in ''.join(codons) or '?' in ''.join(codons):
                removed_bases += 3
                continue

            # remove stop codon
            if i == len(sequences[1][0]) - 3:
                removed_bases += 3
                continue

            trimmed_sequences = [''.join(x) for x in zip(trimmed_sequences, codons)]

        # skip sequences that have been 100% reduced
        if len(trimmed_sequences[0]) == 0:
            continue

        with open(out_phy, 'w') as new_phy:
            print('\t{}\t{}'.format(align_data[0], len(trimmed_sequences[0])), file=new_phy)

            for i in range(0, len(sequences[0])):
                spp = sequences[0][i]
                spp_trim_seq = trimmed_sequences[i]
                print(spp, file=new_phy)
                for block in range(0, len(spp_trim_seq), 60):
                    print(spp_trim_seq[block: block + 60], file=new_phy)

if __name__ == '__main__':
    main()
