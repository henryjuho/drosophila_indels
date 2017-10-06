#!/usr/bin/env python

from __future__ import print_function
import anavar_utils as an
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


def reformat_mle(line, n_classes, var_type, converged, b_hit, model, p, dfe):

    """
    converts mle estimate line from anavar results file to multiple lines in long form
    :param line: dict
    :param n_classes: int
    :param var_type: str
    :param converged: bool
    :param b_hit: list
    :param model: str
    :param p: int
    :param dfe: str
    :return: list
    """

    # determine variant string
    if var_type == 'indel':
        variants = ['ins_', 'del_']
    else:
        variants = ['']

    # prepare header
    out_lines = []
    new_header = ['run', 'imp', 'exit_code',
                  'theta', 'scale', 'shape', 'gamma', 'e',
                  'var_type', 'site_class', 'sel_type', 'lnL',
                  'boundaries_hit', 'converged', 'model', 'params', 'AIC']

    non_variable_vals = {'run': line['run'], 'imp': line['imp'], 'exit_code': line['exit_code'],
                         'lnL': line['lnL'], 'boundaries_hit': '|'.join(b_hit), 'converged': converged,
                         'model': model, 'params': p, 'AIC': aic(p, line['lnL'])}

    # process neutral variants then selected variants
    for sel in ['neu_', 'sel_']:

        # only ever one class of neutral variants
        if sel == 'neu_':
            n_class = 1
        # if sel then number of classes given
        else:
            n_class = n_classes

        for site_class in range(0+1, n_class+1):

            # ie for ins and del of just snp
            for variant in variants:
                if variant == '':
                    var_str = 'snp'
                else:
                    var_str = variant.rstrip('_')

                # get value for each column in output table
                row = []
                for col in new_header:

                    # these values same for each row
                    if col in non_variable_vals.keys():
                        col_val = non_variable_vals[col]

                    elif col == 'var_type':
                        col_val = var_str

                    elif col == 'site_class':
                        col_val = site_class

                    elif col == 'sel_type':
                        col_val = sel.rstrip('_')

                    # leaving 'theta', 'scale', 'shape', 'gamma', 'e'
                    else:
                        if dfe == 'continuous' and sel == 'sel_':
                            class_str = ''
                        else:
                            class_str = '_' + str(n_class)

                        key_str = '{}{}{}{}'.format(sel, variant, col, class_str)

                        try:
                            col_val = line[key_str]
                        except KeyError:
                            col_val = 'na'

                    row.append(col_val)

                out_lines.append(row)

    return new_header, out_lines


def main():

    counter = 0

    all_res = []
    for res in sys.stdin:

        res = res.rstrip()

        if 'equal_t' in res:
            constraint = '_equal_t'
        else:
            constraint = ''

        results = an.ResultsFile(open(res))
        mle = results.ml_estimate()
        n_class = results.num_class()
        variant = results.data_type()
        dfe = results.dfe()
        free_params = len(results.free_parameters())

        mod_name = '{}_{}_{}class{}'.format(variant, dfe, n_class, constraint)

        reformed = reformat_mle(mle, n_class, variant, results.converged(),
                                results.bounds_hit(gamma_r=(-1e5, 1e3), theta_r=(1e-10, 0.1)),
                                mod_name, free_params, dfe)

        if counter == 0:
            all_res.append(reformed[0])
            counter += 1

        for x in reformed[1]:
            all_res.append(x)

    # calc AIC
    aics = sorted([z[-1] for z in all_res[1:]])
    best_aic = aics[0]

    print(*all_res[0] + ['delta_AIC'], sep=',')

    for line in all_res[1:]:
        delta = delta_aic(line[-1], best_aic)
        print(*line + [delta], sep=',')

if __name__ == '__main__':
    main()
