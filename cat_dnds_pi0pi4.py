#!/usr/bin/env python

from __future__ import print_function
import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d', help='File with dnds estimates', required=True)
parser.add_argument('-p', help='File with pi0 pi4 estimates', required=True)
args = parser.parse_args()

# variables
pi_vals = {x.split()[0]: [x.split()[1], x.split()[2]] for x in open(args.p) if not x.startswith('trans_id')}

print('trans', 'length', 'dN', 'dS', 'pi0', 'pi4', sep='\t')
for line in open(args.d):
    if not line.startswith('gene'):
        trans_id = line.split()[0]
        try:
            trans_pi = pi_vals[trans_id]
        except KeyError:
            continue
        out = line.rstrip() + '\t' + '\t'.join(trans_pi)
        print(out)
