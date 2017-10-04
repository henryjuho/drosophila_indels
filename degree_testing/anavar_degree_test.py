#!/usr/bin/env python

from qsub import q_sub

vcf = ('/fastdata/bop15hjb/drosophila_data/dmel/analysis_ready_data/'
       'dmel_17flys.gatk.raw.snps.exsnpindel.recalibrated.filtered_'
       't95.0.pass.dpfiltered.50bp_max.bial.rmarked.polarised.annotated.ar.degen.vcf.gz')

callsites = '/fastdata/bop15hjb/drosophila_data/dmel_ref/dmel.callablesites.summary_with_degen.csv'

out_dir = '/fastdata/bop15hjb/drosophila_data/dmel/anavar/degree_variation/'

for degree in [25, 50, 75, 100, 150, 200, 300]:

    out = '{}dmel_cds_v_4fold_snps_continuous_dfe_degree{}'.format(out_dir, degree)
    cmd = ('cds_vs_neutral_anavar_snps.py '
           '-vcf {} '
           '-n 17 -c 1 -dfe continuous '
           '-call_csv {}  '
           '-neu_type 4fold '
           '-out_pre {} -degree {}'
           '').format(vcf, callsites, out, degree)

    q_sub([cmd], out=out, t=48)
