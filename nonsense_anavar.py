#!/usr/bin/env python

from __future__ import print_function
import anavar_utils as an
import argparse
from qsub import q_sub
from vcf2raw_sfs import vcf2sfs
from cds_vs_neutral_anavar import read_callable_csv, sfs2counts
from nonsense_stats import prem_freqs_call, gather_chromo_prems


def prepare_nonsense_snp_sfs(vcf, call, n, sel_sfs, sel_m):

    """
    gets sfs from vcf and prepares as anavar input
    :param vcf: str
    :param call: dict
    :param n: int
    :param sel_sfs: list
    :param sel_m: int
    :return: dict
    """

    # extract site frequencies
    neu_sfs = vcf2sfs(vcf_name=vcf, mode='snp',
                      auto_only=True, skip_hetero=True,
                      degen=4)

    # convert to correct format for anavar
    sfs_sel = sfs2counts(sel_sfs, n)
    sfs_neu = sfs2counts(neu_sfs, n)

    # get callable sites
    neu_m = call['ALL']['fourfold']['pol']

    # construct control file sfs
    sfs_m = {'selected_SNP': (sfs_sel, sel_m), 'neutral_SNP': (sfs_neu, neu_m)}

    return sfs_m


def sel_v_neu_anavar_nonsense(vcf, call, constraint, n, c, dfe, alg, nnoimp, maximp,
                              out_stem, search, degree, spread, evolgen, prem_files):

    """
    submits anavar jobs to cluster after writing required files etc
    :param vcf: str
    :param call: dict
    :param constraint: str
    :param n: int
    :param c: int
    :param dfe: str
    :param alg: str
    :param nnoimp: int
    :param maximp: int
    :param out_stem: str
    :param search: int
    :param degree: int
    :param spread: int
    :param evolgen: bool
    :param prem_files: list
    :return: None
    """

    anavar_path = '/shared/evolgen1/shared_data/program_files/sharc/'

    anavar_cmd = '{path}anavar1.22 {ctl} {rslts} {log} {seed}'

    # sort file names
    ctl_name = out_stem + '.control.txt'

    # get nonsense data in
    nonsense_dict = gather_chromo_prems(prem_files)
    sel_sfs, sel_m = prem_freqs_call(nonsense_dict)

    # make control file
    sfs_data = prepare_nonsense_snp_sfs(vcf, call, n, sel_sfs, sel_m)
    ctl = an.SNPNeuSelControlFile()

    ctl.set_alg_opts(search=search, alg=alg, key=3,
                     epsabs=1e-20, epsrel=1e-9, rftol=1e-9,
                     maxtime=3600, optional=True,
                     maximp=maximp, nnoimp=nnoimp)

    ctl.set_data(sfs_data, n, dfe=dfe, c=c, gamma_r=(-5e4, 1e3), theta_r=(1e-10, 0.1), r_r=(0.01, 100),
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

            #  seed = random.randint(1, 1e6)
            seed = i

            split_stem = '{}.split{}'.format(out_stem, i)

            result_name = split_stem + '.results.txt'
            log_name = split_stem + '.log.txt'

            print(result_name, file=res_list)

            # call anavar
            rep_cmd = anavar_cmd.format(path=anavar_path, ctl=ctl_name, rslts=result_name, log=log_name, seed=seed)

            q_sub([rep_cmd], out=split_stem, jid=split_stem.split('/')[-1] + '.sh', t=8, evolgen=evolgen)
            hjids.append(split_stem.split('/')[-1] + '.sh')

    # hold job to merge outputs
    merge_out = out_stem + '.merged.results.txt'
    gather = 'cat {} | gather_searches.py {}'.format(res_file_list, merge_out)
    q_sub([gather], out=out_stem + '.merge', hold=hjids, evolgen=evolgen)


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-nonsense_list', help='list file of nonsense data files per chromo', required=True)
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
                        choices=['NLOPT_LD_SLSQP', 'NLOPT_LD_LBFGS', 'NLOPT_LN_NELDERMEAD'])
    parser.add_argument('-nnoimp', help='nnoimp value', default=1, type=int)
    parser.add_argument('-maximp', help='maximp value', default=3, type=int)
    parser.add_argument('-split', help='Number of jobs to split runs across, each job will run the control file once'
                                       'with a different seed given to anavar', default=1, type=int)
    parser.add_argument('-degree', help='changes degree setting in anavar', default=50, type=int)
    parser.add_argument('-out_pre', help='File path and prefix for output', required=True)
    parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
    args = parser.parse_args()

    # variables
    call_site_dict = read_callable_csv(args.call_csv)
    out_pre = args.out_pre
    nonsense_files = [x.rstrip() for x in open(args.nonsense_list)]

    # construct process
    sel_v_neu_anavar_nonsense(vcf=args.vcf, call=call_site_dict,
                              constraint=args.constraint,
                              n=args.n, c=args.c, dfe=args.dfe,
                              alg=args.alg,
                              nnoimp=args.nnoimp, maximp=args.maximp,
                              out_stem=out_pre,
                              search=args.n_search,
                              degree=args.degree,
                              spread=args.split,
                              evolgen=args.evolgen,
                              prem_files=nonsense_files)

if __name__ == '__main__':
    main()
