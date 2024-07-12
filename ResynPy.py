#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on '7/03/19'

@author: 'kalexiou'
"""

import argparse
import os
import sys
import os.path as op
import datetime
from plot_pair_genotypes import *
import cProfile
import pstats
from compare_pairs import *


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Analysis of F2 genotyping data for detection of highly complementary individuals.")

    parser.add_argument('--genos', metavar='TAB FILE', required=True,
                        help='Tab-delimited file (.tab) containing the '
                             'genotyping data of the F2 population. Markers '
                             'should be in columns and individuals in '
                             'rows. Genotypes should be in the format '
                             '"A,B,H,-", where "-" represents missing data. This argument has to be used together '
                             'with --markers. [Required]')
    parser.add_argument('--markers', metavar='TAB FILE', required=True,
                        help='A 2-column tab-delimited file with the markers used for the F2 genotyping, in the format '
                             'of "chromosome<tab>marker name". Marker order per chromosome should match that of '
                             'markers in the genotyping file. [Required]')
    parser.add_argument('--scores_file', metavar='FILE', required=True,
                        help='A tab-delimited file containing user-defined scores for '
                             'the different combinations of genotypes, assigned '
                             'during the comparison of the individuals. If not one available, the user can use the '
                             'scores_default.tab file found in the downloaded github directory. [Required]')
    parser.add_argument('--results_dir', metavar='STR', default='ResynPy_results',
                        help='Name of the results directory [Default: ResynPy_results]')
    parser.add_argument('--not_phased', action='store_true', help='Use this argument if your genotyping data are not '
                                                                  'phased. [Default: FALSE]')
    parser.add_argument('--invariable', metavar='FLOAT', default=0.01,
                        help='Keep individual pairs that have a ratio of AA or '
                             "BB combinations that are lower than FLOAT. [Default: 0.01]")
    parser.add_argument('--score_threshold', metavar='FLOAT', default=0.8,
                        help='Keep individual pairs with total score ratio higher than FLOAT. [Default: 0.8]')
    parser.add_argument('--hetero', metavar='FLOAT', default=0.5,
                        help="Keep individuals that have heterozygosity ratio lower than FLOAT. "
                             '[Default: 0.5]')
    parser.add_argument('--out_prefix', metavar='STR', help='Prefix for the output files.[Default: resynpyOut]',
                        default='resynpyOut')
    parser.add_argument('--cprofile', metavar='STR', help='')

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
        phased = 'no'
    else:
        phased = 'yes'

    run_analysis(param.genos,
                 param.markers,
                 param.results_dir,
                 scoresdict(param.scores_file),
                 param.hetero,
                 param.invariable,
                 param.score_threshold,
                 param.out_prefix,
                 phased)

    return


def plot_pairs(param):

    pairs = '/'.join((param.results_dir, param.out_prefix+'_selected-pairs.tab'))
    fig_pref = param.results_dir+'/'+param.out_prefix+'_top10pairs'
    
    ind_dict = indv2genos(param.genos)
    generate_graph(ind_dict, pairs, fig_pref, param.markers)

    logfile(param, '------\n')
    logfile(param, 'Step 2: Figure generation for the 10 pairs of individuals with the highest score\n')
    logfile(param, '------\n\n')

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


def scoresdict(s_file):
    scores_dict = od()

    with open(s_file) as f:
        for line in f:
            geno, score = line.rstrip('\n\r').split('\t')
            scores_dict[''.join(sorted(geno))] = float(score)

    return scores_dict


def main():
    arguments, variable = parse_arguments()

    if not op.exists('%s' % arguments.results_dir):
        os.makedirs('%s' % arguments.results_dir)

    if (arguments.genos and not arguments.markers) or (arguments.markers and not arguments.genos):
        sys.exit('\nExiting process...\nArguments --markers and --genos should be used together\n')

    pair_comparison(arguments, variable)

    with open(op.join(arguments.results_dir, arguments.out_prefix+'_selected-pairs.tab')) as f:
        f1 = f.readlines()
    
    if len(f1) != 1:
        if sys.platform == 'win32':
            cmd = ' '.join(('python (or python3) plot_pair_genotypes.py', '/'.join((arguments.results_dir,
                                                                                    arguments.out_prefix +
                                                                                    '_selected-pairs.tab')),
                            arguments.genos, arguments.markers,
                            arguments.results_dir + '/' + arguments.out_prefix + '_top10pairs'))

            print('\nYou are using a win32 operating system!\n\n')
            print('\tPlease remember to order manually file {}+_selected-pairs.tab, in order to have the best '
                  'pairs on the top of the list.\n\n'.format(arguments.out_prefix))
            print('\tAfter ordering {}+_selected-pairs.tab file, run the following command to generate the graphs for '
                  'the top10 pairs:\n\n'.format(arguments.out_prefix))
            print('\t\t %s' % cmd)
        
            end = datetime.datetime.now()
            logfile(arguments, '--------\n')
            logfile(arguments, 'End time: ' + str(end) + '\n')
            logfile(arguments, '--------\n\n')

        else:
            plot_pairs(arguments)
            
            end = datetime.datetime.now()
            logfile(arguments, '--------\n')
            logfile(arguments, 'End time: ' + str(end) + '\n')
            logfile(arguments, '--------\n\n')
    else:
        print('\n\t!!! Sorry...No pairs could be selected with the parameters used. Consider allowing higher '
              'percentage of heterozygosity in your individuals and/or higher percentage of invariable sites. !!!\n')
        logfile(arguments, "\t!!! Sorry...No pairs could be selected with the parameters used. Consider allowing "
                           "higher percentage of heterozygosity in your individuals and/or higher percentage of "
                           "invariable sites. !!!\n\n")
        
        end = datetime.datetime.now()
        logfile(arguments, '--------\n')
        logfile(arguments, 'End time: ' + str(end) + '\n')
        logfile(arguments, '--------\n\n')
        

if __name__ == '__main__':
    argus, __ = parse_arguments()
    profiler = cProfile.Profile()
    profiler.enable()
    main()
    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats('cumtime')
    stats.dump_stats(argus.cprofile)
    stats.strip_dirs()
