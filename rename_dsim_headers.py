#!/usr/bin/env python

from __future__ import print_function
import sys

for line in sys.stdin:
    if not line.startswith('>'):
        print(line.rstrip('\n'))
    else:
        if line.startswith('>NODE'):
            print('_'.join(line.split('_')[0:2]))
        else:
            print(line.rstrip('\n'))
