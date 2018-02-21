#!/usr/bin/env python

from __future__ import print_function
from summary_stats import pi, theta_w, tajimas_d
import pysam
import argparse
from cds_vs_neutral_anavar import sfs2counts
import subprocess


def bed_call_sites(call_fa, region_bed):

    """
    returns number of callable sites for specied region in window
    :param call_fa: pysam.FastaFile()
    :param region_bed: pysam.TabixFile()
    :return: int
    """

    contigs = region_bed.contigs
    call_sites = 0

    for contig in contigs:
        call_seq = call_fa.fetch(contig)
        regions = region_bed.fetch(contig, parser=pysam.asTuple())

        for reg in regions:
            call_sub_seq = call_seq[int(reg[1]): int(reg[2])]
            call_sites += call_sub_seq.count('K')

    return call_sites


def main():

    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf')
    parser.add_argument('-fourfold', help='bed file of fourfold sites with final column listing gene name',
                        required=True)
    parser.add_argument('-call_fa')
    args = parser.parse_args()

    # files
    vcf = args.vcf
    bed = pysam.TabixFile(args.fourfold)
    call = pysam.FastaFile(args.call_fa)

    # get callable
    n_call = bed_call_sites(call_fa=call, region_bed=bed)

    print('callable: {}'.format(n_call))

    # get site freqs
    freq_cmd = ('bedtools intersect -header -a {} -b {} | '
                '~/sfs_utils/vcf2raw_sfs.py -mode snp -auto_only -skip_hetero'
                '').format(vcf, args.fourfold)

    freq_list = subprocess.Popen(freq_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    freq_list = [float(x) for x in freq_list]

    str_freq = sfs2counts(freq_list, 17)
    print('sfs: {}'.format(str_freq))

    # get summary stats
    pi_val = pi(17, freq_list) / float(n_call)
    tw = theta_w(17, len(freq_list)) / float(n_call)
    tajd = tajimas_d(17, freq_list)

    print('pi: {}\ntheta_w: {}\nTajD: {}'.format(pi_val, tw, tajd))


if __name__ == '__main__':
    main()
