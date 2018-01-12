#!/usr/bin/env python

from __future__ import print_function
import anavar_utils as an
import argparse
from qsub import q_sub
import random
from collections import Counter
from vcf2raw_sfs import vcf2sfs

# import gzip
# import pysam


# def cds_coords(cds_bed, auto_only):
#
#     """
#     takes a bed file of cds coords in the format chromo\tstart\tstop\ttrans1,trans2,trans3\n
#     :param cds_bed: str
#     :param auto_only: bool
#     :return: dict
#     """
#
#     gene_cds = {}
#     sex_chr = {'X', 'XHet', 'Y', 'YHet'}
#
#     for line in gzip.open(cds_bed):
#         chromo, start, stop, gene = line.split()
#         gene = gene.split(',')[0].split('-')[0]
#
#         if gene not in gene_cds:
#             # exclude sex chr if specified
#             if auto_only and chromo in sex_chr:
#                 continue
#
#             gene_cds[gene] = [chromo, []]
#
#         gene_cds[gene][1].append((int(start), int(stop)))
#
#     return gene_cds
#
#
# def cds_sfs(gene_bed, call_fa, vcf, mode):
#
#     """
#
#     :param gene_bed:
#     :param call_fa:
#     :param vcf:
#     :param mode:
#     :return:
#     """
#
#     all_genes = cds_coords(gene_bed, auto_only=True)
#     for gene in all_genes.keys():
#
#         gene_data = all_genes[gene]
#         gene_chromo = gene_data[0]
#         coords = gene_data[1]
#         chromo_seq = pysam.FastaFile(call_fa).fetch(gene_chromo)
#
#         gene_freqs = []
#         gene_callable = 0
#
#         for coord in coords:
#             # get callable_sites
#             n_callable = chromo_seq[coord[0]: coord[1]].count('K')
#             gene_callable += n_callable
#
#             # get sfs
#             freqs = vcf2sfs(vcf_name=vcf,
#                             mode=mode,
#                             chromo=gene_chromo, start=coord[0], stop=coord[1],
#                             skip_hetero=True)
#             gene_freqs += list(freqs)
#
#         yield gene_freqs, gene_callable
#
#
# def neutral_sfs():
#     pass


def read_callable_csv(csv):

    """
    reads in a callable sites summary file
    :param csv: str
    :return: dict
    """

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

    """
    resamples sfs with replacement
    :param site_freqs: list
    :return: list
    """

    resamp_sfs = []
    for i in range(0, len(site_freqs)):
        random_no = random.randint(0, len(site_freqs) - 1)
        resamp_sfs.append(site_freqs[random_no])
    return resamp_sfs


def sfs2counts(freq_list, n):

    """
    converts list of sfs to a condensed sfs that can be plotted or passed to anavar etc
    :param freq_list: list
    :param n: int
    :return: list
    """

    pos_biallelic_freqs = [round(i / float(n), 3) for i in range(1, int(n))]

    counted_sorted_sfs = sorted(Counter([str(x) for x in freq_list]).most_common(), key=lambda z: z[0])
    sfs_freq_dict = {x[0]: x[1] for x in counted_sorted_sfs}

    counts = []
    for frequency in pos_biallelic_freqs:
        try:
            no_var = sfs_freq_dict[str(frequency)]
        except KeyError:
            no_var = 0  # account for freqs with 0 variants

        counts.append(no_var)

    return counts


def prepare_indel_sfs(vcf, call, n):

    """
    gets sfs from vcf and prepares as anavar input
    :param vcf: str
    :param call: dict
    :param n: int
    :return: dict
    """

    # extract site frequencies
    del_sfs = vcf2sfs(vcf_name=vcf, mode='del',
                      auto_only=True, skip_hetero=True,
                      regions=['CDS_frameshift', 'CDS_non_frameshift'])

    n_d_sfs = vcf2sfs(vcf_name=vcf, mode='del',
                      auto_only=True, skip_hetero=True,
                      regions=['intergenic', 'intron'])

    ins_sfs = vcf2sfs(vcf_name=vcf, mode='ins',
                      auto_only=True, skip_hetero=True,
                      regions=['CDS_frameshift', 'CDS_non_frameshift'])

    n_i_sfs = vcf2sfs(vcf_name=vcf, mode='ins',
                      auto_only=True, skip_hetero=True,
                      regions=['intergenic', 'intron'])

    # convert to correct format for anavar
    sfs_i = sfs2counts(ins_sfs, n)
    sfs_d = sfs2counts(del_sfs, n)
    sfs_ni = sfs2counts(n_i_sfs, n)
    sfs_nd = sfs2counts(n_d_sfs, n)

    # get callable sites
    sel_m = call['ALL']['CDS']['pol']
    neu_m = call['ALL']['intergenic']['pol'] + call['ALL']['intron']['pol']

    # construct control file sfs
    sfs_m = {'selected_INS': (sfs_i, sel_m), 'selected_DEL': (sfs_d, sel_m),
             'neutral_INS': (sfs_ni, neu_m), 'neutral_DEL': (sfs_nd, neu_m)}

    return sfs_m


def prepare_snp_sfs(vcf, call, n):

    """
    gets sfs from vcf and prepares as anavar input
    :param vcf: str
    :param call: dict
    :param n: int
    :return: dict
    """

    # extract site frequencies
    sel_sfs = vcf2sfs(vcf_name=vcf, mode='snp',
                      auto_only=True, skip_hetero=True,
                      regions=['CDS_frameshift', 'CDS_non_frameshift'])

    neu_sfs = vcf2sfs(vcf_name=vcf, mode='snp',
                      auto_only=True, skip_hetero=True,
                      degen=4)

    # convert to correct format for anavar
    sfs_sel = sfs2counts(sel_sfs, n)
    sfs_neu = sfs2counts(neu_sfs, n)

    # get callable sites
    sel_m = call['ALL']['CDS']['pol']
    neu_m = call['ALL']['fourfold']['pol']

    # construct control file sfs
    sfs_m = {'selected_SNP': (sfs_sel, sel_m), 'neutral_SNP': (sfs_neu, neu_m)}

    return sfs_m


def sel_v_neu_anavar(mode, vcf, call, constraint, n, c, dfe, alg, out_stem, search, degree, spread, evolgen):

    """
    submits anavar jobs to cluster after writing required files etc
    :param mode: str
    :param vcf: str
    :param call: dict
    :param constraint: str
    :param n: int
    :param c: int
    :param dfe: str
    :param alg: str
    :param out_stem: str
    :param search: int
    :param degree: int
    :param spread: int
    :param evolgen: bool
    :return: None
    """

    anavar_path = '/shared/evolgen1/shared_data/program_files/sharc/'

    anavar_cmd = '{path}anavar1.22 {ctl} {rslts} {log} {seed}'

    # sort file names
    ctl_name = out_stem + '.control.txt'

    # make control file
    if mode == 'snp':
        sfs_data = prepare_snp_sfs(vcf, call, n)
        ctl = an.SNPNeuSelControlFile()

    else:
        sfs_data = prepare_indel_sfs(vcf, call, n)
        ctl = an.IndelNeuSelControlFile()

    ctl.set_alg_opts(search=search, alg=alg, key=3,
                     epsabs=1e-20, epsrel=1e-9, rftol=1e-9,
                     maxtime=3600, optional=True)

    ctl.set_data(sfs_data, n, dfe=dfe, c=c, gamma_r=(-5e4, 1e2), theta_r=(1e-10, 0.1), r_r=(0.01, 100),
                 scale_r=(0.1, 5000.0))
    if degree != 50:
        ctl.set_dfe_optional_opts(degree=degree, optional=True)
    ctl.set_constraint(constraint)
    ctl_contents = ctl.construct()
    with open(ctl_name, 'w') as control:
        control.write(ctl_contents)

    res_file_list = out_stem + '.allres.list.txt'
    hjids = []
    with open(res_file_list, 'w') as res_list:

        # split into requested jobs
        for i in range(1, spread+1):

            split_stem = '{}.split{}'.format(out_stem, i)

            result_name = split_stem + '.results.txt'
            log_name = split_stem + '.log.txt'

            print(result_name, file=res_list)

            # call anavar
            rep_cmd = anavar_cmd.format(path=anavar_path, ctl=ctl_name, rslts=result_name, log=log_name, seed=i)

            q_sub([rep_cmd], out=split_stem, jid=split_stem.split('/')[-1] + '.sh', t=8, evolgen=evolgen)
            hjids.append(split_stem.split('/')[-1] + '.sh')

    # hold job to merge outputs
    merge_out = out_stem + '.merged.results.txt'
    gather = 'cat {} | gather_searches.py {}'.format(res_file_list, merge_out)
    q_sub([gather], out=out_stem + '.merge', hold=hjids, evolgen=evolgen)


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', help='variant type to run on', choices=['snp', 'indel'], required=True)
    parser.add_argument('-call_csv', help='Callable sites summary file', required=True)
    parser.add_argument('-vcf', help='VCF file to extract site frequencies from', required=True)
    parser.add_argument('-n', help='Sample size', required=True)
    parser.add_argument('-c', help='Number of classes to run model with', required=True, type=int)
    parser.add_argument('-dfe', help='type of dfe to fit, discrete or continuous', default='discrete',
                        choices=['discrete', 'continuous'])
    parser.add_argument('-constraint', help='Constraint for model', choices=['none', 'equal_mutation_rate'],
                        default='none')
    parser.add_argument('-n_search', help='Number of searches to conduct per job', default=500, type=int)
    parser.add_argument('-alg', help='Algorithm to use', default='NLOPT_LD_SLSQP',
                        choices=['NLOPT_LD_SLSQP', 'NLOPT_LD_LBFGS'])
    parser.add_argument('-split', help='Number of jobs to split runs across, each job will run the control file once'
                                       'with a different seed given to anavar', default=1, type=int)
    parser.add_argument('-degree', help='changes degree setting in anavar', default=50, type=int)
    parser.add_argument('-out_pre', help='File path and prefix for output', required=True)
    parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
    args = parser.parse_args()

    # variables
    call_site_dict = read_callable_csv(args.call_csv)
    out_pre = args.out_pre

    # construct process
    sel_v_neu_anavar(mode=args.mode,
                     vcf=args.vcf, call=call_site_dict,
                     constraint=args.constraint,
                     n=args.n, c=args.c, dfe=args.dfe,
                     alg=args.alg,
                     out_stem=out_pre,
                     search=args.n_search,
                     degree=args.degree,
                     spread=args.split,
                     evolgen=args.evolgen)

if __name__ == '__main__':
    main()
