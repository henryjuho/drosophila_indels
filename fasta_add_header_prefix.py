#!/usr/bin/env python

import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-fa', help='Fasta file to add prefix to headers', required=True)
parser.add_argument('-pre', help='Prefix', required=True)
parser.add_argument('-truncate', help='If specified will truncate header at first whitespace', action='store_true',
                    default=False)
args = parser.parse_args()

# variables
fa = args.fa
prefix = args.pre
trunc = args.truncate
out = fa.rstrip('.fa') + '.rename.fa'

# rename
with open(out, 'w') as new_fa:
    for line in open(fa):
        if line.startswith('>'):
            line = line.replace('>', '>' + prefix)
            if trunc is True:
                line = line.split()[0] + '\n'
            new_fa.write(line)
        else:
            new_fa.write(line)

print('Done')
