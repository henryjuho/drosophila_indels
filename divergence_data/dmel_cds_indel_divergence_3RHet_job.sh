#!/bin/bash

source ~/.bash_profile

#$-l h_rt=8:00:00
#$-l mem=6G
#$-l rmem=2G

#$-P evolgen
#$-q evolgen.q

#$-N dmel_cds_indel_divergence_3RHet_job.sh
#$-o /home/bop15hjb/drosophila_indels/summary_stats/dmel_cds_indel_divergence_3RHet.out
#$-e /home/bop15hjb/drosophila_indels/summary_stats/dmel_cds_indel_divergence_3RHet.error

#$-V

./indel_divergence.py -wga /fastdata/bop15hjb/drosophila_data/wga/multiple_alignment/dmel.dsim.dyak.wga.bed.gz -bed /fastdata/bop15hjb/drosophila_data/dmel_ref/dmel_cds.bed -out /home/bop15hjb/drosophila_indels/summary_stats/dmel_cds_indel_divergence.txt -chromo 3RHet
