from __future__ import print_function
import argparse
import os
from qsub import q_sub


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-cds_fa_dir', help='cds fasta directory', required=True)
    parser.add_argument('-out_dir', help='output directory', required=True)
    parser.add_argument('-vcf', help='SNP vcf path', required=True)
    parser.add_argument('-call_fa', help='Callable sites fasta file', required=True)
    parser.add_argument('-evolgen', help='if specified will run on lab queue', default=False, action='store_true')
    args = parser.parse_args()

    out_files = []
    out_dir = args.out_dir
    autos = ('2R', '2RHet', '2L', '2LHet', '3R', '3RHet', '3L', '3LHet', '4')

    # per chromo jobs for extracting nonsense data
    for x in [args.cds_fa_dir + y for y in os.listdir(args.cds_fa_dir)]:
        if not x.endswith('.fasta.gz'):
            continue

        chromo = x.split('-')[1]

        if chromo not in autos:
            continue

        outstem = x.split('/')[-1].replace('.fasta.gz', '')
        out = out_dir + outstem + '.premstops.txt'

        out_files.append(out)

        extract_cmd = ('./extract_prem_stops.py '
                       '-cds_fa {cds_fa} '
                       '-chr {chromo} '
                       '-vcf {vcf} '
                       '-call_fa {c_fa} '
                       '-n 17 '
                       '-unfolded '
                       '-out {output}').format(
            cds_fa=x,
            chromo=chromo,
            vcf=args.vcf,
            c_fa=args.call_fa,
            output=out)

        q_sub([extract_cmd], out=out_dir + outstem, evolgen=args.evolgen)

    # write list file
    list_file = out_dir + 'chromo_nonsense_list.txt'

    with open(list_file) as list_out:
        print(*out_files, sep='\n', file=list_out)


if __name__ == '__main__':
    main()
