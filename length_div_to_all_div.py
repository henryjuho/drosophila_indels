from __future__ import print_function
import sys

data_dict = {'deletion': 0, 'insertion': 0}

lines = [x.rstrip() for x in sys.stdin]
calls = 0

for line in lines:
    line = line.rstrip()
    if line.startswith('length'):
        continue
    else:
        n_var, var_type = int(line.split('\t')[1]), line.split('\t')[-1]
        data_dict[var_type] += n_var
        calls = int(line.split('\t')[2])

print('<51', data_dict['deletion'], calls, data_dict['deletion']/float(calls), 'deletion', sep='\t')
print('<51', data_dict['insertion'], calls, data_dict['insertion']/float(calls), 'insertion', sep='\t')
