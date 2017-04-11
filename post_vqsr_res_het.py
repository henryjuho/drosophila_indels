#!/usr/bin/env python

import sys
from qsub import *

for vcf in sys.stdin:
    dirname = vcf[:vcf.rfind('/')+1]
    out = dirname + vcf.replace('.pass.vcf', '').split('_')[-1] + 'resid_het.txt'
    cmd_line = 'cat ' + vcf + ' | python per_indiv_res_het.py > ' + out
    q_sub([cmd_line], out=out.replace('.txt', ''))
