#!/usr/bin/env python

from __future__ import print_function
import sys

# loop
counter = 0
for line in sys.stdin:
    counter += 1
    if counter <= 3:
        continue
    else:
        line = line.split()
        chromo = line[4]
        start = str(int(line[5]) - 1)
        end = line[6]
        print('\t'.join([chromo, start, end]))
