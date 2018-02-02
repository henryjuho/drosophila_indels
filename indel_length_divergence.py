#!/usr/bin/env python

from __future__ import print_function
import argparse
import subprocess
import sys


def popen_grab(cmd):
    output_lines_list = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    return output_lines_list


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-wga', help='Whole genome alignment bed file', required=True)
    parser.add_argument('-bed', help='Coordinates to calc divergence for, bed format', required=True)
    parser.add_argument('-out', help='Output path and file', required=True)
    args = parser.parse_args()

    with open(args.out, 'w') as out_file:
        print('length\tn_variants\tcallable\tdivergence\tindel_type', file=out_file)

    callable_cmd = 'bedtools intersect -a {} -b {} | ~/WGAbed/wga_bed_summary.py -callable'.format(args.wga, args.bed)

    print(callable_cmd, file=sys.stdout)

    call_sites = popen_grab(callable_cmd)[0].split('\t')

    n_sites = int(call_sites[1])

    # per length per indel type
    for var_type in ['deletion', 'insertion']:
        for length in ['1', '2', '3', '4+', '6+']:

            if length == '4+':
                lens = range(4, 51, 3) + range(5, 51, 3)
                lens = ','.join(sorted([str(x) for x in lens]))
            elif length == '6+':
                lens = ','.join(sorted([str(x) for x in range(6, 51, 3)]))
            else:
                lens = length

            indel_cmd = ('bedtools intersect -a {} -b {} | '
                         '~/WGAbed/wga_bed_indels.py '
                         '-max_length 50 -min_coverage 3 '
                         '-ref_specific -lengths {} | '
                         '~/WGAbed/polarise_wga_ref_indels.py '
                         '-indel_type {} | grep -v ^X | grep -v ^Y | wc -l').format(
                args.wga, args.bed, lens, var_type)

            print(indel_cmd, file=sys.stdout)

            n_indels = int(popen_grab(indel_cmd)[0])

            if n_sites == 0:
                div = 0.0
            else:
                div = float(n_indels)/float(n_sites)

            with open(args.out, 'a') as out_file:
                print(length, n_indels, n_sites, div, var_type, file=out_file, sep='\t')


if __name__ == '__main__':
    main()
