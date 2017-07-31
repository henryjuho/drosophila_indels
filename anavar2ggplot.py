#!/usr/bin/env python

import sys


def reformat_line(line, header):

    out_lines = []
    new_header = ['run', 'imp', 'exit_code', 'theta', 'gamma', 'e', 'var_type', 'site_class', 'lnL', 'rep', 'region']
    header = header.split()
    out_data = {x: '' for x in new_header}
    data_to_sort = {}
    for i in range(0, len(header)):
        current_col = header[i]
        current_val = line.split()[i]
        if 'ins' not in current_col and 'del' not in current_col:
            out_data[current_col] = current_val
        else:
            data_to_sort[current_col] = current_val

    # detrimine if class 1 or 2 is largest by ins gamma
    ga = data_to_sort['ins_gamma_1']
    gb = data_to_sort['ins_gamma_2']

    if float(ga) > float(gb):
        c1 = '1'
        c2 = '2'
    else:
        c1 = '2'
        c2 = '1'

    # write multiple rows in long form
    for var_type in ['ins', 'del']:
        for site_class in [(c1, '1'), (c2, '2')]:
            out_data['theta'] = data_to_sort['{}_theta_{}'.format(var_type, site_class[0])]
            out_data['gamma'] = data_to_sort['{}_gamma_{}'.format(var_type, site_class[0])]
            out_data['e'] = data_to_sort['{}_e_{}'.format(var_type, site_class[0])]
            out_data['var_type'] = var_type
            out_data['site_class'] = site_class[1]

            out_lines.append('\t'.join([out_data[y] for y in new_header]))

    return '\t'.join(new_header), out_lines


def main():

    header_line = ''
    data_line = 0
    for line in sys.stdin:
        if line.startswith('run'):
            header_line = line
        else:
            data_line += 1
            new_lines = reformat_line(line, header_line)

            if data_line == 1:
                print new_lines[0]

            print '\n'.join(new_lines[1])

if __name__ == '__main__':
    main()
