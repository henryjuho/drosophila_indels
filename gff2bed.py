#!/usr/bin/env python

from __future__ import print_function
import sys

for line in sys.stdin:
    if line.startswith('#'):
        continue
    else:
        line = line.split()
        chromo, start, end = line[0], int(line[3])-1, line[4]

        print(chromo, start, end, sep='\t')
