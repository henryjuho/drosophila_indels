#!/usr/bin/env python

from __future__ import print_function
import os

directory = '/fastdata/bop15hjb/drosophila_data/dmel/anavar/degree_variation'

no = 0
for x in [directory + f for f in os.listdir(directory) if f.endswith('.allreps.results.txt')]:
    res_contents = open(x).readlines()
    degree = x.split('degree')[-1].split('.')[0]
    if no == 0:
        print(','.join(res_contents[0].split() + ['degree']))
        no += 1
    print(','.join(res_contents[1].split() + [degree]))
