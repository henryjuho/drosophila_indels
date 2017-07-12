#!/usr/bin/env python

from __future__ import print_function
import argparse
import random
from qsub import *


def no_trans_longest(transcript):
    transcript = sorted(transcript)
    long_len = transcript[-1][0]

    counter = 0
    for x in transcript[::-1]:
        if x[0] == long_len:
            counter += 1
        else:
            break

    return counter


def control_file(seq_file, out_file, tree_file_in, ctl_file):

    control_str = ('      seqfile = {}\n'
                   '     treefile = {}\n'
                   '      outfile = {}           * main result file name\n'
                   '\n'
                   '        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen\n'
                   '      verbose = 1  * 0: concise; 1: detailed, 2: too much\n'
                   '      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic\n'
                   '                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise\n'
                   '\n'
                   '      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs\n'
                   '    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n'
                   '\n'
                   '*        ndata = 10\n'
                   '        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n'
                   '       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\n'
                   '*   aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)\n'
                   '                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own\n'
                   '\n'
                   '        model = 0\n'
                   '                   * models for codons:\n'
                   '                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches\n'
                   '                   * models for AAs or codon-translated AAs:\n'
                   '                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F\n'
                   '                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)\n'
                   '\n'
                   '      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;\n'
                   '                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;\n'
                   '                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;\n'
                   '                   * 13:3normal>0\n'
                   '\n'
                   '        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below\n'
                   '        Mgene = 0\n'
                   '                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff\n'
                   '                   * AA: 0:rates, 1:separate\n'
                   '\n'
                   '    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated\n'
                   '        kappa = 2  * initial or fixed kappa\n'
                   '    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate\n'
                   '        omega = .4 * initial or fixed omega, for codons or codon-based AAs\n'
                   '\n'
                   '    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha\n'
                   '        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)\n'
                   '       Malpha = 0  * different alphas for genes\n'
                   '        ncatG = 8  * # of categories in dG of NSsites models\n'
                   '\n'
                   '        getSE = 0  * 0: dont want them, 1: want S.E.s of estimates\n'
                   ' RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)\n'
                   '\n'
                   '   Small_Diff = .5e-6\n'
                   '    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?\n'
                   '*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed\n'
                   '       method = 0  * Optimization method 0: simultaneous; 1: one branch a time\n'
                   '\n'
                   '* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,\n'
                   '* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt.,\n'
                   '* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt.,\n'
                   '* 10: blepharisma nu.\n'
                   '* These codes correspond to transl_table 1 to 11 of GENEBANK.\n'
                   '\n').format(seq_file, tree_file_in, out_file)

    with open(ctl_file, 'w') as control:
        print(control_str, file=control)


def tree_file(tree_file_name):

    with open(tree_file_name, 'w') as tree:
        tree_str = ('\t3\t1\n\n'
                    '((dmel, dsim), dyak);')
        print(tree_str, file=tree)


def main():
    # arguments
    parser = argparse.ArgumentParser(description='Runs codeml for each phylip file in the given directory')
    parser.add_argument('-phy_dir',
                        help='Directory of phylip files name in form of gene-transcript.len.phy',
                        required=True)
    parser.add_argument('-out_dir',
                        help='Output directory to write codeml outputs to',
                        required=True)
    args = parser.parse_args()

    # make transcript dict
    phylips = [args.phy_dir + x for x in os.listdir(args.phy_dir) if x.endswith('.phy')]

    tran_dict = {}
    for trans in phylips:
        gene = '-'.join(trans.split('-')[:-1])
        trans_id = trans.split('-')[-1].split('.')[0]
        length = int(trans.split('-')[-1].split('.')[1])
        if gene not in tran_dict.keys():
            tran_dict[gene] = []
        tran_dict[gene].append((length, trans_id, trans))

    # construct codeml job
    cmd_list = []

    counter = 0
    for fly_gene in tran_dict.keys():
        gene_transcripts = sorted(tran_dict[fly_gene])
        no_of_tran_longest = no_trans_longest(gene_transcripts)
        if no_of_tran_longest == 1:
            longest_trans = gene_transcripts[-1]

        # if multiple longest transcripts
        else:
            longest_trans_list = gene_transcripts[-no_of_tran_longest:]

            random_trans_index = random.randint(0, len(longest_trans_list)-1)
            longest_trans = longest_trans_list[random_trans_index]

        longest_phylip = longest_trans[2]

        # prepare codeml files
        basename = args.out_dir + longest_phylip.replace('.phy', '').split('/')[-1]
        tr_file = basename + '.tree.txt'
        tree_file(tr_file)
        codeml_out = basename + '.codeml.out'
        ctrl_file = basename + '.ctl'
        control_file(longest_phylip, codeml_out, tr_file, ctrl_file)

        # designed for parallel on evolgen
        if counter == 20:
            counter = 0
            cmd_list += ['', 'wait', '']

        counter += 1
        codeml_cmd = '/shared/evolgen1/shared_data/program_files/iceberg/codeml ' + ctrl_file + ' &'
        cmd_list.append(codeml_cmd)

    q_sub(cmd_list, out=args.out_dir + 'codeml_runs', t=168, tr=20, evolgen=True)


if __name__ == '__main__':
    main()
