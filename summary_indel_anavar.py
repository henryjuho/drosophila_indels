#!/usr/bin/env python

from __future__ import print_function
import anavar_utils as an
import argparse
from qsub import *
import sys
import random
from multiprocessing import Pool
from collections import Counter


def read_callable_csv(csv):
    call_sites = {}
    for line in open(csv):
        if not line.startswith('contig'):
            info = line.rstrip().split(',')
            contig, reg, all_call, pol_call = info
            if contig not in call_sites.keys():
                call_sites[contig] = {reg: {'all': float(all_call), 'pol': float(pol_call)}}
            else:
                call_sites[contig][reg] = {'all': float(all_call), 'pol': float(pol_call)}
    return call_sites


def resample_replace(site_freqs):
    resamp_sfs = []
    for i in range(0, len(site_freqs)):
        random_no = random.randint(0, len(site_freqs) - 1)
        resamp_sfs.append(site_freqs[random_no])
    return resamp_sfs


def sfs2counts(freq_list, n):
    pos_biallelic_freqs = [round(i / float(n), 3) for i in range(1, int(n))]

    counted_sorted_sfs = sorted(Counter(freq_list).most_common(), key=lambda z: z[0])
    sfs_freq_dict = {x[0]: x[1] for x in counted_sorted_sfs}

    counts = []
    for frequency in pos_biallelic_freqs:
        try:
            no_var = sfs_freq_dict[str(frequency)]
        except KeyError:
            no_var = 0  # account for freqs with 0 variants

        counts.append(no_var)

    return counts


def indel_anavar(arg_tuple):
    ins_sfs, ins_m, del_sfs, del_m, bootstrap, n, region, out_stem = arg_tuple

    anavar_path = '/shared/evolgen1/shared_data/program_files/iceberg/'

    anavar_cmd = '{path}anavar1.2 {ctl} {rslts} {log}'

    results = []

    for i in [0] + range(1, bootstrap+1):

        # sort sfs
        if i == 0:
            sfs_i = ins_sfs
            sfs_d = del_sfs
        else:
            sfs_i = resample_replace(ins_sfs)
            sfs_d = resample_replace(del_sfs)

        # convert to correct format for anavar
        sfs_i = sfs2counts(sfs_i, n)
        sfs_d = sfs2counts(sfs_d, n)

        # sort file names
        ctl_name = out_stem + '.rep{}.control.txt'.format(i)
        result_name = out_stem + '.rep{}.results.txt'.format(i)
        log_name = out_stem + '.rep{}.log.txt'.format(i)

        # construct control file
        sfs_m = {'INS': (sfs_i, ins_m), 'DEL': (sfs_d, del_m)}
        ctl = an.Indel1ControlFile()
        ctl.set_data(sfs_m, n, dfe='discrete', c=2)
        ctl_contents = ctl.construct()
        with open(ctl_name, 'w') as control:
            control.write(ctl_contents)

        # call anavar
        rep_cmd = anavar_cmd.format(path=anavar_path, ctl=ctl_name, rslts=result_name, log=log_name)
        subprocess.call(rep_cmd, shell=True)

        # process results
        with open(result_name) as rep_results:
            results_data = an.ResultsFile(rep_results)
            header = list(results_data.header()) + ['rep', 'region']
            ml_est = results_data.ml_estimate(as_string=True) + '\t{}\t{}'.format(i, region)
            if i == 0:
                results.append('\t'.join(header))
            results.append(ml_est)

    return results


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='VCF file to extract site frequencies from', required=True)
    parser.add_argument('-n', help='Sample size', required=True)
    parser.add_argument('-call_csv', help='Callable sites summary file', required=True)
    parser.add_argument('-bootstrap', help='Number of bootstrap replicates', default=0, type=int)
    parser.add_argument('-out_pre', help='File path and prefix for output', required=True)
    parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
    parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
    args = parser.parse_args()

    # submission loop
    if args.sub is True:
        command_line = [' '.join([x for x in sys.argv if x != '-sub' and x != '-evolgen'])]
        q_sub(command_line, args.out_pre, evolgen=args.evolgen, tr=4, t=96)
        sys.exit(0)

    # variables
    call_site_dict = read_callable_csv(args.call_csv)
    out_pre = args.out_pre
    processes_args = []

    # set up each region run
    for region in ['ALL', 'CDS', 'intergenic', 'intron']:

        file_stem = out_pre + '.' + region

        # extract site frequencies
        del_sfs_cmd = ('/home/bop15hjb/sfs_utils/vcf2raw_sfs.py '
                       '-vcf {vcf} '
                       '-mode del -auto_only -skip_hetero'
                       ).format(vcf=args.vcf)

        ins_sfs_cmd = ('/home/bop15hjb/sfs_utils/vcf2raw_sfs.py '
                       '-vcf {vcf} '
                       '-mode ins -auto_only -skip_hetero'
                       ).format(vcf=args.vcf)

        if region == 'CDS':
            region_flags = ' -region CDS_frameshift -region CDS_non_frameshift'
        elif region != 'ALL':
            region_flags = ' -region {}'.format(region)
        else:
            region_flags = ''

        del_sfs_cmd += region_flags
        ins_sfs_cmd += region_flags

        del_sfs = subprocess.Popen(del_sfs_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
        ins_sfs = subprocess.Popen(ins_sfs_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]

        # get callable sites
        del_m = call_site_dict['ALL'][region]['pol']
        ins_m = call_site_dict['ALL'][region]['pol']

        # construct process
        anavar_args = (ins_sfs, ins_m, del_sfs, del_m, args.bootstrap, args.n, region, file_stem)
        processes_args.append(anavar_args)

    # initiate pool and run processes
    anavar_pool = Pool(processes=len(processes_args))
    region_results = anavar_pool.map(indel_anavar, processes_args)

    # main out
    with open(out_pre + '.allregions.results.txt', 'w') as main_out:
        # write results
        counter = 0
        for reg in region_results:

            counter += 1

            # printers header once
            if counter == 1:
                print(reg[0], file=main_out)

            # trims off header
            reg_data = reg[1:]

            for line in reg_data:
                print(line, file=main_out)


if __name__ == '__main__':
    main()
