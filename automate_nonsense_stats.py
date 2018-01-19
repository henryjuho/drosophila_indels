#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import q_sub
import os


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-cda_fa_dir', help='cds fasta directory', required=True)
    parser.add_argument('-out_dir', help='output directory', required=True)
    parser.add_argument('-vcf', help='SNP vcf path', required=True)
    parser.add_argument('-call_fa', help='Callable sites fasta file', required=True)
    parser.add_argument('-evolgen', help='if specified will run on lab queue', default=False, action='store_true')
    args = parser.parse_args()

    out_files = []
    jobs = []
    autos = ('2R', '2RHet', '2L', '2LHet', '3R', '3RHet', '3L', '3LHet', '4')

    # per chromo jobs
    for x in [args.cds_fa_dir + y for y in os.listdir(args.cds_fa_dir)]:

        chromo = x.split('-')[1]

        if chromo not in autos:
            continue

        outstem = x.split('/')[-1].replace('.fasta.gz', '')
        out = args.out_dir + outstem + '.premstops.txt'
        job = outstem + '.sh'

        out_files.append(out)
        jobs.append(job)

        extract_cmd = ('./extract_prem_stops.py '
                       '-cds_fa {cds_fa} '
                       '-chr {chromo} '
                       '-vcf {vcf} '
                       '-call_fa {c_fa} '
                       '-n 17 '
                       '-out {output}').format(
            cds_fa=x,
            chromo=chromo,
            vcf=args.vcf,
            c_fa=args.call_fa,
            output=out)

        q_sub([extract_cmd], out=args.out_dir + outstem, jid=job, evolgen=args.evolgen)

    # gather outputs and calc stats
    gather_cmd = ('./nonsense_stats.py ' +
                  ' '.join(['-infile ' + x for x in out_files]) +
                  ' -out_file ' + args.out_dir + 'dmel_nonsense_stats.txt')

    q_sub([gather_cmd], out=args.out_dir + 'dmel_nonsense_stats', hold=jobs, evolgen=args.evolgen)


if __name__ == '__main__':
    main()
