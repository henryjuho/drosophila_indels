#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-dir', help='Directory containing maf files to ensure single coverage for', required=True)
parser.add_argument('-ref_name', help='Name of reference species', required=True)
parser.add_argument('-evolgen', help='If specified will submit to lab queue', default=False, action='store_true')
args = parser.parse_args()

# variables
directory = args.dir
maf_list = [maf for maf in os.listdir(directory) if maf.endswith('.maf')]
out_dir = directory + 'single_coverage/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
ref = args.ref_name
evolgen = args.evolgen

# run single_cov2
for maf in maf_list:
    in_maf = directory + maf
    output = out_dir + maf.rstrip('.maf') + '.sing.maf'
    cmd_line = ('single_cov2 ' +
                in_maf + ' [R=' + ref + '] > ' +
                output)
    q_sub([cmd_line], out=out_dir + 'single_cov2_' + maf.rstrip('.maf'),
          jid='single_cov2_' + maf.rstrip('.maf') + '.sh', evolgen=evolgen)
