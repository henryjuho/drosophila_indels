#!/usr/bin/env python

from __future__ import print_function
import sys


def aic(n_param, max_lnl):

    """
    calculates AIC as desicribed here: https://en.wikipedia.org/wiki/Akaike_information_criterion
    :param n_param: int
    :param max_lnl: float
    :return:
    """

    k = n_param
    lnl = max_lnl

    return (2 * k) - (2 * lnl)


def main():
    # r params + neutral params + sel indel params
    no_params = {'sel_v_neu_continuous': 15 + 3 + 8,
                 'sel_v_neu_2class': 15 + 3 + 12,
                 'sel_v_neu_1class': 15 + 3 + 6,
                 'sel_v_neu_continuous_equal_t': 15 + 3 + 7,
                 'sel_v_neu_2class_equal_t': 15 + 3 + 11,
                 'sel_v_neu_1class_equal_t': 15 + 3 + 5}

    for line in sys.stdin:
        if line.startswith('run'):
            print(line.rstrip(), 'n_params', 'AIC', sep=',')
        else:
            split_line = line.rstrip().split(',')
            model = split_line[14]
            max_l = float(split_line[11])
            n = no_params[model]
            line_aic = aic(n, max_l)

            print(line.rstrip(), n, line_aic, sep=',')

if __name__ == '__main__':
    main()
