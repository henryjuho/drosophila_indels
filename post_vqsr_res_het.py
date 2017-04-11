#!/usr/bin/env python

import sys
from qsub import *

for vcf in sys.stdin:
    vcf = vcf.rstrip('\n')
    dirname = vcf[:vcf.rfind('/')+1]
    out = dirname + vcf.replace('.pass.vcf', '').split('_')[-1] + '_resid_het.txt'
    cmd_line = 'cat ' + vcf + ' | python per_indiv_res_het.py > ' + out
    q_sub([cmd_line], out=out.replace('.txt', ''))
