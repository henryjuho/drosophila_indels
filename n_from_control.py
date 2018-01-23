from __future__ import print_function
import sys

print('length\tn_segsites\tvar_type')
for line in sys.stdin:

    n = sum([int(x) for x in line.split(':')[-1].split(',')])
    var_type = line.split(':')[-2].split('_')[1]
    length = line.split(':')[-3].split('_len')[-1].split('.')[0]

    print(length, n, var_type, sep='\t')
