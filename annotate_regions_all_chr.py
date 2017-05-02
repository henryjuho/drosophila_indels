#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gff', help='GFF file to read annotations from, note squence names must match VCF', required=True)
parser.add_argument('-vcf', help='VCF file to annotate variants in', required=True)
parser.add_argument('--feat_merge_gene_proxy', help='If specified will merge all features to get gene coords,'
                                                    'only valid if only gene features annotated in gff',
                    action='store_true', default=False)
parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
args = parser.parse_args()

# variables
gff = args.gff
vcf = args.vcf
gene_prox = args.feat_merge_gene_proxy
evolgen = args.evolgen

# get chromosome list
grep_cmd = 'grep -v ^# ' + vcf + ' | cut -f 1 | uniq'
chromo_list = subprocess.Popen(grep_cmd, stdout=subprocess.PIPE, shell=True).communicate()[0].split('\n')[:-1]
chromo_list = [x for x in chromo_list if not x.startswith('NODE')]

# loop through chromo list and submit annotation job for each
vcf_outs = []
hold_list = []
for chromo in chromo_list:
    chr_jid = 'annotation_' + chromo + '.sh'
    hold_list.append(chr_jid)
    annotate_vcf_cmd = ('vcf_region_annotater.py '
                        '-gff ' + gff + ' '
                        '-vcf ' + vcf + ' '
                        '-chr ' + chromo + ' ')
    if gene_prox is True:
        annotate_vcf_cmd += '--feat_merge_gene_proxy'
    out_vcf = vcf.replace('.vcf', '.annotated.' + chromo + '.vcf')
    vcf_outs.append(out_vcf)
    q_sub([annotate_vcf_cmd],
          out=out_vcf.rstrip('.vcf'),
          jid=chr_jid,
          evolgen=evolgen)

# submit concat job
cat_cmd = ('catVCFs.py '
           '-out_vcf ' + vcf.replace('.vcf', '.annotated.vcf'))
for out_vcf in vcf_outs:
    cat_cmd += ' -vcf ' + out_vcf

q_sub([cat_cmd],
      out=vcf.replace('.vcf', '.annotated'),
      hold=hold_list,
      evolgen=evolgen)
