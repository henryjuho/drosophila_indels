#!/bin/bash

source ~/.bash_profile

#$-l h_rt=8:00:00
#$-l mem=6G
#$-l rmem=2G

#$-P evolgen
#$-q evolgen.q

#$-N dmel_del_ar.sfs_job.sh
#$-o indel_sfs/dmel_del_ar.sfs.out
#$-e indel_sfs/dmel_del_ar.sfs.error

#$-V

~/sfs_utils/vcf2raw_sfs.py -vcf /fastdata/bop15hjb/drosophila_data/dmel/analysis_ready_data/dmel_17flys.gatk.raw.indels.recalibrated.filtered_t95.0.pass.dpfiltered.50bp_max.bial.rmarked.polarised.annotated.ar.vcf.gz -mode del -auto_only -skip_hetero -region intergenic_ar -region intron_ar | sort | uniq -c | while read i; do echo " "$i; done | tr -s " " > indel_sfs/dmel_del_ar.sfs.txt
