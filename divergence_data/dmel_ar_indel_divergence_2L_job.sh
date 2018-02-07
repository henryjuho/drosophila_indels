#!/bin/bash

source ~/.bash_profile

#$-l h_rt=8:00:00
#$-l mem=6G
#$-l rmem=2G

#$-P evolgen
#$-q evolgen.q

#$-N dmel_ar_indel_divergence_2L_job.sh
#$-o /home/bop15hjb/drosophila_indels/divergence_data/dmel_ar_indel_divergence_2L.out
#$-e /home/bop15hjb/drosophila_indels/divergence_data/dmel_ar_indel_divergence_2L.error

#$-V

./indel_divergence.py -wga /fastdata/bop15hjb/drosophila_data/wga/multiple_alignment/dmel.dsim.dyak.wga.bed.gz -bed /fastdata/bop15hjb/drosophila_data/wga/multiple_alignment/dmel_ancestral_repeats.wga.bed.gz -out /home/bop15hjb/drosophila_indels/divergence_data/dmel_ar_indel_divergence.txt -chromo 2L
