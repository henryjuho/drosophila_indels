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


def delta_aic(current_d, best_d):

    """
    returns delta_AIC, which is defined as the difference between the AIC for the best fitting model
    and the model on the row of the spreadsheet - from Kai's email 14th Sept 17
    :param current_d: float
    :param best_d: float
    :return:
    """

    return best_d - current_d


def main():
    # r params + neutral params + sel indel params
    no_params = {'sel_v_neu_continuous': {'snp': 15 + 2 + 4, 'ins': 15 + 4 + 8, 'del': 15 + 4 + 8},
                 'sel_v_neu_3class': {'snp': 15 + 2 + 9, 'ins': 15 + 4 + 18, 'del': 15 + 4 + 18},
                 'sel_v_neu_2class': {'snp': 15 + 2 + 6, 'ins': 15 + 4 + 12, 'del': 15 + 4 + 12},
                 'sel_v_neu_1class': {'snp': 15 + 2 + 3, 'ins': 15 + 4 + 6, 'del': 15 + 4 + 6},
                 'sel_v_neu_continuous_equal_t': {'snp': 15 + 1 + 4, 'ins': 15 + 2 + 8, 'del': 15 + 2 + 8},
                 'sel_v_neu_3class_equal_t': {'snp': 15 + 1 + 9, 'ins': 15 + 2 + 18, 'del': 15 + 2 + 18},
                 'sel_v_neu_2class_equal_t': {'snp': 15 + 1 + 6, 'ins': 15 + 2 + 12, 'del': 15 + 2 + 12},
                 'sel_v_neu_1class_equal_t': {'snp': 15 + 1 + 3, 'ins': 15 + 2 + 6, 'del': 15 + 2 + 6}}

    contents = [line for line in sys.stdin]

    print(contents[0].rstrip(), 'n_params', 'AIC', 'deltaAIC', sep=',')

    results = []
    aics = []
    for res in contents[1:]:
        split_line = res.rstrip().split(',')
        model = split_line[14]
        var_type = split_line[8]
        max_l = float(split_line[11])
        n = no_params[model][var_type]
        line_aic = aic(n, max_l)
        split_line += [n, line_aic]
        results.append(split_line)
        aics.append(line_aic)

    best_aic = sorted(aics)[0]

    for line in results:

        d_aic = delta_aic(line[-1], best_aic)

        print(','.join([str(x) for x in line]), d_aic, sep=',')

if __name__ == '__main__':
    main()
