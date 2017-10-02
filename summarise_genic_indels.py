#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import *
import sys
import gzip
import pysam
from summary_stats import pi, theta_w, tajimas_d


def bed_to_dict(c_bed):

    gene_dict = {}

    for line in gzip.open(c_bed):
        chromo, start, stop, gene_id = line.rstrip().split('\t')
        gene_id = gene_id.split(',')

        for trans_id in gene_id:
            if chromo not in gene_dict.keys():
                gene_dict[chromo] = {}
            if trans_id not in gene_dict[chromo].keys():
                gene_dict[chromo][trans_id] = []

            gene_dict[chromo][trans_id] += range(int(start), int(stop))

    return gene_dict


def main():
    # arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-call_fa',
                        help='Callable sites fasta file',
                        required=True)
    parser.add_argument('-vcf',
                        help='Vcf file to extract site frequencies from',
                        required=True)
    parser.add_argument('-cds_bed',
                        help='Bed file of zerofold sites, in form chr\tstart\tstop\gene_id',
                        required=True)
    parser.add_argument('-out',
                        help='Output file location and name',
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
        q_sub(command_line,
              out=args.out.replace('.txt', '') + 'gene_pi0_pi4',
              evolgen=args.evolgen, t=48,
              mem=15, rmem=15)
        sys.exit()

    # variables
    call_fa = pysam.FastaFile(args.call_fa)
    vcf = pysam.VariantFile(args.vcf)
    gene_coords = bed_to_dict(args.cds_bed)
    number_samples = len(vcf.header.samples)
    out = open(args.out, 'w')

    # gene by gene calcs
    print('trans_id', 'pi_indel', 'theta_indel', 'tajd_indel', sep='\t', file=out)

    for chromosome in gene_coords.keys():
        chr_string = call_fa.fetch(chromosome)
        for trans in gene_coords[chromosome].keys():

            call_sites = ''
            allele_freqs = []
            for pos in gene_coords[chromosome][trans]:

                # get callable site
                call_pos = chr_string[pos]
                call_sites += call_pos

                # get vcf site (try to)
                var_record = [x for x in vcf.fetch(chromosome, pos, pos+1)]
                if len(var_record) == 1:
                    allele_freq = round(var_record[0].info['AC'][0] / float(number_samples * 2), 3)
                    allele_freqs.append(allele_freq)

            # count callable sites for transcript
            n_callable = call_sites.upper().count('K')

            # calc pi
            pie = pi(number_samples, allele_freqs)
            theta = theta_w(number_samples, len(allele_freqs))
            tajd = tajimas_d(number_samples, allele_freqs)

            if n_callable != 0:
                pie_per_site = pie / float(n_callable)
                theta_per_site = theta / float(n_callable)
                tajd_per_site = tajd / float(n_callable)
            else:
                pie_per_site, theta_per_site, tajd_per_site = 0.0, 0.0, 0.0

            print(trans, pie_per_site, theta_per_site, tajd_per_site, sep='\t', file=out)

    out.close()

if __name__ == '__main__':
    main()
