#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on '7/03/19'

@author: 'kalexiou'
"""

import argparse
import os
import os.path as op
import datetime
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Analysis of F2 genotyping data for detection of highly complementary individuals.")

    parser.add_argument('--genos', metavar='TAB FILE',
                        help='Tab-delimited file (.tab) containing the '
                             'genotyping data of the F2 population. Markers '
                             'should be in columns and individuals in '
                             'rows. Genotypes should be in the format '
                             '"A,B,H,-", where "-" represents missing data. This argument has to be used together '
                             'with --markers')

    parser.add_argument('--markers', metavar='TAB FILE',
                        help='A 2-column tab-delimited file with the markers used for the F2 genotyping, in the format '
                             'of "chromosome<tab>marker name". Marker order per chromosome should match that of '
                             'markers in the genotyping file.')

    parser.add_argument('--results_dir', metavar='STR', default='ResynPy_results',
                        help='Name of the results directory [Default: ResynPy_results]')

    parser.add_argument('--not_phased', action='store_true', help='Use this argument if your genotyping data are not '
                                                                  'phased. [Default: FALSE]')

    parser.add_argument('--scores_file', metavar='FILE', default='scores_default.tab',
                        help='A tab delimited file containing user-defined scores for '
                             'the different combinations of genotypes, assigned '
                             'during the comparison of the individuals. Nucleotides in'
                             'genotypes should be separated with a "/". '
                             'e.g.: A/H<tab>0.75. [Default: scores_default.tab]')

    parser.add_argument('--invariable', metavar='FLOAT', default=0.1,
                        help='Keep individual pairs that have a percentage of AA or '
                             "BB combinations that is lower than the argument's value. [Default: 0.1]")

    parser.add_argument('--hetero', metavar='FLOAT', default=0.5,
                        help="Keep individuals that have a heterozygosity lower than the argument's value. "
                             '[Default: 0.5]')

    parser.add_argument('--out_prefix', metavar='STR', help='Prefix for the output files.[Default: resynpyOut]',
                        default='resynpyOut')

    variables = vars(parser.parse_args())

    return parser.parse_args(), variables


def logfile(param, string):
    with open('%s/PairComparison.log' % param.results_dir, 'a') as log:
        log.writelines(string)

    return


def pair_comparison(param, arguments):

    now2 = datetime.datetime.now()

    logfile(param, '==============================\n')
    logfile(param, 'Analysis of F2 genotyping data\n')
    logfile(param, '==============================\n\n')

    logfile(param, '---------\n')
    logfile(param, 'variables:\n')
    logfile(param, '---------\n\n')

    for k, v in arguments.items():
        logfile(param, '\t' + k + ': ' + str(v) + '\n')

    logfile(param, '\n----------\n')
    logfile(param, 'Start time: ' + str(now2) + '\n')
    logfile(param, '----------\n\n')
    logfile(param, '------\n')
    logfile(param, 'Step 1: Pairwise comparison of individuals\n')
    logfile(param, '------\n\n')
    if param.not_phased:
        cmd = ' '.join(('python3 compare_pairs.py', param.genos, param.markers, param.results_dir,
                        param.scores_file, str(param.hetero), str(param.invariable), param.out_prefix, 'no'))
    else:
        cmd = ' '.join(('python3 compare_pairs.py', param.genos, param.markers, param.results_dir,
                        param.scores_file, str(param.hetero), str(param.invariable), param.out_prefix, 'yes'))

    logfile(param, '\t' + cmd + '\n\n')
    os.system(cmd)

    return


def plot_pairs(param):

    cmd = ' '.join(('python3 plot_pair_genotypes.py', '/'.join((param.results_dir,
                                                                param.out_prefix+'_selected-pairs.tab')),
                    param.genos, param.markers, param.results_dir+'/'+param.out_prefix+'_top10pairs'))

    os.system(cmd)

    logfile(param, '------\n')
    logfile(param, 'Step 2: Figure generation for the 10 pairs of individuals with the highest score\n')
    logfile(param, '------\n\n')

    logfile(param, '\t' + cmd + '\n\n')

    end = datetime.datetime.now()

    logfile(param, '--------\n')
    logfile(param, 'End time: ' + str(end) + '\n')
    logfile(param, '--------\n\n')

    return


def roh_write_variables(param, args):
    logfile(param, '---------\n')
    logfile(param, 'variables:\n')
    logfile(param, '---------\n\n')

    for k, v in args.items():
        to_none = ['scores_file', 'invariable', 'hetero', 'fig_prefix']
        if k in to_none:
            v = 'None'
        else:
            v = v

        logfile(param, '\t' + k + ': ' + str(v) + '\n')


def main():
    arguments, variable = parse_arguments()

    if not op.exists('%s' % arguments.results_dir):
        os.makedirs('%s' % arguments.results_dir)

    if (arguments.genos and not arguments.markers) or (arguments.markers and not arguments.genos):
        sys.exit('\nExiting process...\nArguments --markers and --genos should be used together\n')

    pair_comparison(arguments, variable)
    plot_pairs(arguments)


if __name__ == '__main__':
    main()
