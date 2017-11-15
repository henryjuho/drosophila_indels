from __future__ import print_function
import sys

del_dat = {x.split()[0]: x.split()[1:] for x in open(sys.argv[1]) if not x.startswith('x')}
ins_dat = {x.split()[0]: x.split()[1:] for x in open(sys.argv[2]) if not x.startswith('x')}

print('x', 'pN', 'pS', sep='\t')
for x in sorted(del_dat.keys()):
    print(x, int(del_dat[x][0]) + int(ins_dat[x][0]), int(del_dat[x][1]) + int(ins_dat[x][1]), sep='\t')
