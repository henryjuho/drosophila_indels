from __future__ import print_function
import sys

cds_dat = {x.split()[1]: x.split()[0] for x in open(sys.argv[1])}
noncode_dat = {x.split()[1]: x.split()[0] for x in open(sys.argv[2])}

print('x', 'pN', 'pS', sep='\t')
for x in sorted(cds_dat.keys()):
    print(x, cds_dat[x], noncode_dat[x], sep='\t')
