#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: 'kalexiou'
"""

import argparse
import itertools
import os
import sys
import os.path as op
import datetime
import multiprocessing as mp
from plot_pair_genotypes import *
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
    parser.add_argument('--scores_file', metavar='FILE', default='ResynPy/scores_default.tab',
                        help='A tab-delimited file containing user-defined scores for '
                             'the different combinations of genotypes, assigned '
                             'during the comparison of the individuals. If not one available, the user can use the '
                             'scores_default.tab file found in the downloaded github directory. '
                             '[Default: ResynPy/scores_default.tab]')
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
    parser.add_argument('--num_cpus', metavar='INT', type=int, default=4,
                        help='Number of CPUs to use for the parallelization. [Default: 4]')

    variables = vars(parser.parse_args())

    return parser.parse_args(), variables


def logfile(param, string):
    with open('%s/PairComparison.log' % param.results_dir, 'a') as log:
        log.writelines(string)

    return


def plot_pairs(param):
    pairs = '/'.join((param.results_dir, param.out_prefix + '_selected-pairs.tab'))
    fig_pref = param.results_dir + '/' + param.out_prefix + '_top10pairs'

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


def split_list(lst, chunk_size):
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]


def run_analysis_wrapper(args):
    return run_analysis(*args)


def parallel_analysis(lst, chunk_size, new_content, marker_file, ideal_score, scoresd, invar, total_score_threshold,
                      phased):
    chunks = list(split_list(lst, chunk_size))
    with mp.Pool() as pool:
        results = pool.map(run_analysis_wrapper, [(index_chunk, new_content, marker_file,
                                                   ideal_score, scoresd, invar, total_score_threshold,
                                                   phased) for index_chunk in chunks])
    return results


def main():
    arguments, variables = parse_arguments()

    if not op.exists('%s' % arguments.results_dir):
        os.makedirs('%s' % arguments.results_dir)

    if (arguments.genos and not arguments.markers) or (arguments.markers and not arguments.genos):
        sys.exit('\nExiting process...\nArguments --markers and --genos should be used together\n')

    now2 = datetime.datetime.now()

    logfile(arguments, 'Start time: ' + str(now2) + '\n')
    logfile(arguments, '\n')
    logfile(arguments, '==============================\n')
    logfile(arguments, 'Analysis of F2 genotyping data\n')
    logfile(arguments, '==============================\n\n')

    logfile(arguments, '---------\n')
    logfile(arguments, 'variables:\n')
    logfile(arguments, '---------\n\n')

    for k, v in variables.items():
        logfile(arguments, '\t' + k + ': ' + str(v) + '\n')

    logfile(arguments, '\n')
    logfile(arguments, '------------------------------------------\n')
    logfile(arguments, 'Step 1: Pairwise comparison of individuals\n')
    logfile(arguments, '------------------------------------------\n\n')

    if arguments.not_phased:
        phased = 'no'
    else:
        phased = 'yes'

    with open(arguments.markers) as tmp:
        tmp1 = tmp.readlines()
        ideal_score: int = len(tmp1)

    with open(arguments.genos, "rt") as f:
        f1 = f.readlines()
        content = [x.rstrip('\n\r').split('\t') for x in f1]

        new_content, disc_indv_list, disc_individuals = filter_hetero(content, arguments.hetero)

        logfile(arguments, '\tComparison Statistics' + '\n')
        logfile(arguments, '\t---------------------\n\n')
        logfile(arguments, '\t\ttotal number of markers      : ' + str(ideal_score) + '\n')
        logfile(arguments, '\t\ttotal number of individuals  : ' + str(len(content[1:])) + '\n')
        logfile(arguments, '\t\tindividuals not meeting\n'
                           '\t\t  heterozygosity threshold   : ' + str(disc_individuals) + '\n')
        logfile(arguments, '\t\tindividuals to compare       : ' + str(len(new_content[1:])) + '\n')
        logfile(arguments, '\t\tpairs to compare             : ' +
                str((len(new_content[1:]) * (len(new_content[1:]) - 1)) / 2).split('.')[0] + '\n')

        index_pairs_list = list(itertools.combinations(range(1, len(new_content)), 2))
        final_output = parallel_analysis(lst=index_pairs_list,
                                         chunk_size=int(len(index_pairs_list) / arguments.num_cpus),
                                         new_content=new_content,
                                         marker_file=arguments.markers, ideal_score=ideal_score,
                                         scoresd=scoresdict(arguments.scores_file), invar=arguments.invariable,
                                         total_score_threshold=arguments.score_threshold, phased=phased)
        outfile = op.join(arguments.results_dir, arguments.out_prefix + '_selected-pairs.tab')
        temp_file = op.join(arguments.results_dir, arguments.out_prefix + '.temp')

    with open(temp_file, 'w') as out:
        for i in final_output:
            for z in i[0]:
                out.writelines('\t'.join(z) + '\n')

    a = pd.read_table(temp_file, header=None, sep='\t')

    a.columns = ['#indiv_pair', 'score', 'invariable_pairs', 'hetero1|hetero2', 'similarity_to_hybrid',
                 'recomb_no1', 'recomb_no2', 'sum_recomb']
    b = a.sort_values(by=["score", "sum_recomb"], ascending=[False, True])
    b.to_csv(outfile, index=False, sep='\t')

    os.remove(temp_file)

    accepted_pairs = sum([i[1] for i in final_output])
    discarded_pairs_list = list()
    for i in final_output:
        for z in i[2]:
            discarded_pairs_list.append(z)
    
    discarded_pairs_score = sum([i[3] for i in final_output])
    discarded_pairs_invar = sum([i[4] for i in final_output])

    logfile(arguments, '\t\tdiscarded pairs due to score : ' +
            str(discarded_pairs_score) + '\n')
    logfile(arguments, '\t\tpairs not meeting invariable\n'
                       '\t\t  site threshold             : ' +
            str(discarded_pairs_invar) + '\n')
    logfile(arguments, '\t\tACCEPTED PAIRS               : ' +
            str(accepted_pairs) + '\n\n')

    with open(op.join(arguments.results_dir, arguments.out_prefix + '_discarded_individuals-pairs.tab'), 'w') as disc:
        discarded_data = disc_indv_list + discarded_pairs_list
        for a in discarded_data:
            disc.writelines('\t'.join(a) + '\n')

    # plot top-10 individuals
    with open(outfile) as f:
        f1 = f.readlines()
        if len(f1) != 1:
            print('plotting top pairs')
            plot_pairs(arguments)
            print('plotted pairs')
        else:
            print('\n\t!!! Sorry...No pairs could be selected with the parameters used. Consider allowing higher '
                  'percentage of heterozygosity in your individuals and/or higher percentage of invariable sites. !!!\n')
            logfile(arguments, "\t!!! Sorry...No pairs could be selected with the parameters used. Consider allowing "
                               "higher percentage of heterozygosity in your individuals and/or higher percentage of "
                               "invariable sites. !!!\n\n")

        end = datetime.datetime.now()
        logfile(arguments, 'End time: ' + str(end) + '\n')


if __name__ == '__main__':
    main()
