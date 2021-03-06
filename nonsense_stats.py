#!/usr/bin/env python

from __future__ import print_function
import argparse
import random
from summary_stats import tajimas_d, theta_w, pi


def prem_freqs_call(prem_dict):

    """
    gets allele freq and callable sites form prem dict
    :param prem_dict: dict
    :return: list, int
    """

    call = 0
    minor_freq = []

    for long_trans in prem_dict.keys():

        n_call = int(prem_dict[long_trans][3])
        call += n_call
        freqs = prem_dict[long_trans][4].split(',')
        if len(freqs) == 1:
            continue
        else:
            minor_freq += [float(x) for x in freqs]

    return minor_freq, call


def gather_chromo_prems(in_files):

    """
    merges multiple chromo nonsense data and retains longest trans only for each gene
    :param in_files: list
    :return: dict
    """

    all_data = {}

    for non_file in in_files:

        for line in open(non_file):

            gene, transcript, chromo, trans_len, n_call, allele_freqs = line.rstrip().split('\t')

            if gene not in all_data.keys():
                all_data[gene] = [transcript, chromo, trans_len, n_call, allele_freqs]

            else:
                entered_len = int(all_data[gene][2])
                current_len = int(trans_len)

                if current_len > entered_len:
                    all_data[gene] = [transcript, chromo, trans_len, n_call, allele_freqs]

                elif current_len < entered_len:
                    continue

                else:
                    choice = random.randint(0, 1)

                    if choice == 0:
                        all_data[gene] = [transcript, chromo, trans_len, n_call, allele_freqs]
                    else:
                        continue

    return all_data


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-in_file', help='chromosomal nonsense data', required=True, action='append')
    parser.add_argument('-out_file', help='output file for summary stats', required=True)
    args = parser.parse_args()

    # make dataset
    all_data = gather_chromo_prems(args.in_file)

    # consolidate callable and allele frequencies
    minor_freq, call = prem_freqs_call(all_data)

    # stats
    pi_value = pi(17, minor_freq) / float(call)
    tw = theta_w(17, len(minor_freq)) / float(call)
    tajd = tajimas_d(17, minor_freq)

    with open(args.out_file, 'w') as out:
        print('region\tbin\ttype\tseg_sites\tcallable\ttheta_w\tt_lwr\tt_upr\tpi\t'
              'pi_lwr\tpi_upr\ttajD\ttajD_lwr\ttajD_upr', file=out)
        print('gwide\tnonsense\tsnp\t{seg}\t{call}\t{t}\t0\t0\t{p}\t0\t0\t{taj}\t0\t0'.format(
            seg=len(minor_freq), call=call, t=tw, p=pi_value, taj=tajd), file=out)


if __name__ == '__main__':
    main()
