#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-dir', help='directory containing codeml output files with suffix "codeml.out"', required=True)
args = parser.parse_args()

# grab out files
codeml_outs = [x for x in os.listdir(args.dir) if x.endswith('codeml.out')]

os.chdir(args.dir)

print('gene\tlength\tdN\tdS')
for run in codeml_outs:
    gene = run.split('.')[0]
    length = run.split('.')[1]

    grep_cmd = 'grep "tree length for" {}'.format(run.replace('(', '\(').replace(')', '\)').replace("'", "\\'"))

    dndsdata = subprocess.Popen(grep_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    dn = float(dndsdata[0].split()[-1])
    ds = float(dndsdata[1].split()[-1])

    print(gene, length, dn, ds, sep='\t')
