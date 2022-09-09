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
from collections import OrderedDict as od


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Detection of ROH regions and analysis of F2 genotyping data for detection of highly complementary "
                    "individuals.")

    parser.add_argument('--genos', metavar='TAB FILE',
                        help='Tab-delimited file (.tab) containing the '
                             'genotyping data of the F2 population. Markers '
                             'should be in columns and individuals in '
                             'rows. Genotypes should be in the format '
                             '"A,B,H,-", where "-" represents missing data. Incomapatible with --vcf. This '
                             'argument has to be used together with --markers')

    parser.add_argument('--markers', metavar='TAB FILE',
                        help='A 2-column tab-delimited file with the markers used for the F2 genotyping, in the format '
                             'of "chromosome<tab>marker name". Incombatible with "--vcf"')

    parser.add_argument('--results_dir', metavar='STR', default='ResynPy_results',
                        help='Name of the results directory [Default: ResynPy_results]')

    parser.add_argument('--not_phased', action='store_true', help='Use this argument if your genotyping data are not '
                                                                  'phased. [Default: FALSE]')

    parser.add_argument('--scores_file', metavar='FILE', default='scores_default.tab',
                        help='A tab delimited file containing user-defined scores for '
                             'the different combinations of genotypes, observed '
                             'during the comparison of the individuals. Nucleotides in'
                             'genotypes should be separated with a "/". '
                             'e.g.: A/H<tab>0.75. [Default: scores_default.tab]')

    parser.add_argument('--filter_invariable', metavar='FLOAT', default=0.1,
                        help='Keep individual pairs that have a percentage of AA or '
                             'BB combinations that is lower than the argument value. [Default: 0.1]')

    parser.add_argument('--hetero', metavar='FLOAT', default=0.5,
                        help='Keep individuals that have a heterozygosity lower than the argument value. '
                             '[Default: 0.5]')

    parser.add_argument('--fig_prefix', metavar='STR', default='top10_selected_pairs', help='Prefix to be used for the '
                                                                                            '.png and .pdf figure, '
                                                                                            'showing the genotype '
                                                                                            'profiles of the 10 pairs '
                                                                                            'with the highest score. '
                                                                                            '[Default: top10_selected'
                                                                                            '_pairs]')

    parser.add_argument('--vcf', metavar='VCF FILE', help='Name of the VCF file for the line, for which the ROH '
                                                          'regions will be detected.')

    variables = vars(parser.parse_args())

    return parser.parse_args(), variables


def logfile(param, string):
    with open('%s/bcftools_roh-detection.log' % param.results_dir, 'a') as log:
        log.writelines(string)
    return


def f2logfile(param, string):
    with open('%s/PairComparison.log' % param.results_dir, 'a') as log:
        log.writelines(string)

    return


def roh_bcftools(param):
    com1 = ' '.join(('bcftools roh --AF-dflt 0.4 -G 30 %s' % param.vcf, '>', '%s/ROH_G30' % param.results_dir))
    logfile(param, '\t' + com1 + '\n')
    os.system(com1)

    com2 = ' '.join(('grep -v \'#\' %s/ROH_G30' % param.results_dir, '| awk \'$1=="ST"\'',
                     '|', 'cut -f 3,4,5', '>', '%s/ROH_G30_sites' % param.results_dir))

    logfile(param, '\t' + com2 + '\n')
    os.system(com2)

    with open('sdf', 'w') as out:
        out.writelines('\t'.join(('chromosome', 'position', 'region\n')))
        with open('%s/ROH_G30_sites' % param.results_dir) as f:
            for line in f:
                items = line.rstrip('\n\r').split('\t')
                if items[-1] == '0':
                    region = 'noROH'
                else:
                    region = 'ROH'
                out.writelines('\t'.join((items[0], items[1], region + '\n')))

    com3 = ' '.join(('mv', 'sdf', '%s/ROH_G30_sites' % param.results_dir))
    os.system(com3)
    logfile(param, '\t' + com3 + '\n\n')

    return


def generate_graph(param):
    logfile(param, '------\n')
    logfile(param, 'Step 2: Graph generation\n')
    logfile(param, '------\n\n')
    logfile(param, '\tfigure location: %s/ROH.png\n\n' % param.results_dir)

    x = pd.read_csv('%s/ROH_G30_sites' % param.results_dir, header=0, sep='\t')

    x['position'] = x['position'].div(1000000, axis='index')
    plt.rc('font', size=6)
    g = sns.stripplot(x=x['position'], y=x['chromosome'], jitter=0.2, hue=x['region'], palette="Set2", dodge=True,
                      size=0.6)

    figure = g.get_figure()
    plt.xlabel('position (Mbp)')
    plt.title("Heterozygous (noROH) and homozygous (ROH) regions across the indvidual's genome")
    figure.savefig('%s/ROH.png' % param.results_dir, dpi=400)
    plt.close()

    return


def snps_roh(param):
    logfile(param, '------\n')
    logfile(param, 'Step 3: Generate a file with ROH annotation.\n')
    logfile(param, '------\n\n')

    com = 'head -4 {0}/ROH_G30 | tail -1 | cut -f 3,4,5 > {0}/ROH_G30.header'.format(param.results_dir)
    com1 = ' '.join(('grep -v \'#\' {0}/ROH_G30'.format(param.results_dir), '| awk \'$1=="RG"\'',
                     '| awk \'{OFS="\\t"}{print $3, $4-1, $5}\'',
                     '> {0}/ROH_G30.positions'.format(param.results_dir)))  # I substract 1 nt from the start in
    # order to exclude the start coordinate of the ROH_G30 file.
    com1a = 'cat {0}/ROH_G30.header {0}/ROH_G30.positions > x'.format(param.results_dir)
    com1b = 'mv x {0}/ROH_G30.positions'.format(param.results_dir)

    if param.vcf[-3:] == '.gz':
        com2 = 'vcftools --gzvcf {0} --exclude-bed {1}/ROH_G30.positions --recode ' \
               '--out {1}/noROH_G30'.format(param.vcf, param.results_dir)
        com2a = 'vcftools --gzvcf {0} --bed {1}/ROH_G30.positions --recode ' \
                '--out {1}/ROH_G30'.format(param.vcf, param.results_dir)
    else:
        com2 = 'vcftools --vcf {0} --exclude-bed {1}/ROH_G30.positions --recode ' \
               '--out {1}/noROH_G30'.format(param.vcf, param.results_dir)
        com2a = 'vcftools --vcf {0} --bed {1}/ROH_G30.positions --recode ' \
                '--out {1}/ROH_G30'.format(param.vcf, param.results_dir)

    com3 = ' '.join(('grep -v \'#\' {0}/ROH_G30.recode.vcf'.format(param.results_dir),
                     '| awk \'{OFS="\\t"}{print $1, $2, $4, $5, "ROH"}\'',
                     '> {0}/markers_ROH-anno.tab'.format(param.results_dir)))
    com3a = ' '.join(('grep -v \'#\' {0}/noROH_G30.recode.vcf'.format(param.results_dir),
                      '| awk \'{OFS="\\t"}{print $1, $2, $4, $5, "noROH"}\'',
                      '> {0}/markers_noROH-anno.tab'.format(param.results_dir)))
    com4 = 'echo "chrom\tpos\tref\talt\tregion" > {0}/markers_region-anno.tab'.format(param.results_dir)
    com4a = 'cat {0}/markers_ROH-anno.tab {0}/markers_noROH-anno.tab | sort -k1,1 -k2,2n | uniq >> ' \
            '{0}/markers_region-anno.tab'.format(param.results_dir)

    logfile(param, '\t' + com + '\n')
    logfile(param, '\t' + com1 + '\n')
    logfile(param, '\t' + com1a + '\n')
    logfile(param, '\t' + com1b + '\n')
    logfile(param, '\t' + com2 + '\n')
    logfile(param, '\t' + com2a + '\n')
    logfile(param, '\t' + com3 + '\n')
    logfile(param, '\t' + com3a + '\n')
    logfile(param, '\t' + com4 + '\n')
    logfile(param, '\t' + com4a + '\n')

    for command in [com, com1, com1a, com1b, com2, com2a, com3, com3a, com4, com4a]:
        os.system(command)

    return


def roh_detection(param, arguments):
    now = datetime.datetime.now()

    logfile(param, '===============================================================================================\n')
    logfile(param, 'ROH detection using BCFTOOLS, graph generation and selection of variants inside non-ROH regions\n')
    logfile(param,
            '===============================================================================================\n\n')

    logfile(param, '---------\n')
    logfile(param, 'variables:\n')
    logfile(param, '---------\n\n')

    for k, v in arguments.items():
        logfile(param, '\t' + k + ': ' + str(v) + '\n')

    logfile(param, '\n----------\n')
    logfile(param, 'Start time: ' + str(now) + '\n')
    logfile(param, '----------\n\n')

    logfile(param, '------\n')
    logfile(param, 'Step 1: ROH detection\n')
    logfile(param, '------\n\n')

    roh_bcftools(param)

    generate_graph(param)

    snps_roh(param)

    com = 'rm {0}/markers_ROH-anno.tab {0}/markers_noROH-anno.tab  {0}/ROH_G30.header ' \
          '{0}/ROH_G30_sites {0}/noROH_G30.recode.vcf {0}/noROH_G30.log {0}/ROH_G30.positions ' \
          '{0}/ROH_G30.recode.vcf {0}/ROH_G30.log'.format(param.results_dir)

    logfile(param, '\t' + com + '\n\n')
    os.system(com)

    end = datetime.datetime.now()
    logfile(param, '--------\n')
    logfile(param, 'End time: ' + str(end) + '\n')
    logfile(param, '--------\n\n')

    return


def pair_comparison(param):
    if param.not_phased:
        cmd = ' '.join(('python3 compare_pairs_20220315.py', param.genos, param.markers, param.results_dir,
                        param.scores_file, str(param.hetero), str(param.filter_invariable), 'no'))
    else:
        cmd = ' '.join(('python3 compare_pairs_20220315.py', param.genos, param.markers, param.results_dir,
                        param.scores_file, str(param.hetero), str(param.filter_invariable), 'yes'))

    f2logfile(param, '\t' + cmd + '\n\n')
    os.system(cmd)

    return


def f3_prob(param):
    pass


def f2data_analysis(param, arguments):
    now2 = datetime.datetime.now()

    f2logfile(param, '==============================\n')
    f2logfile(param, 'Analysis of F2 genotyping data\n')
    f2logfile(param, '==============================\n\n')

    f2logfile(param, '---------\n')
    f2logfile(param, 'variables:\n')
    f2logfile(param, '---------\n\n')

    for k, v in arguments.items():
        f2logfile(param, '\t' + k + ': ' + str(v) + '\n')

    f2logfile(param, '\n----------\n')
    f2logfile(param, 'Start time: ' + str(now2) + '\n')
    f2logfile(param, '----------\n\n')
    f2logfile(param, '------\n')
    f2logfile(param, 'Step 1: Pairwise comparison of individuals\n')
    f2logfile(param, '------\n\n')

    pair_comparison(param)

    #    f2logfile(param, '------\n')
    #    f2logfile(param, 'Step 2: Probability calculation for desired combinations in the F3 generation\n')
    #    f2logfile(param, '------\n\n')

    # f3_prob(param)

    return


def plot_pairs(param):

    cmd = ' '.join(('python3 plot_pair_genotypes.py', '/'.join((param.results_dir, 'individual_genetic_distance.tab')),
                    param.genos, param.markers, param.results_dir+'/'+param.fig_prefix))

    os.system(cmd)

    f2logfile(param, '\t' + cmd + '\n\n')

    end = datetime.datetime.now()

    f2logfile(param, '--------\n')
    f2logfile(param, 'End time: ' + str(end) + '\n')
    f2logfile(param, '--------\n\n')
    
    return


def main():
    arguments, variable = parse_arguments()

    if not op.exists('%s' % arguments.results_dir):
        os.makedirs('%s' % arguments.results_dir)

    if arguments.vcf:

        if arguments.vcf[-7:] != '.vcf.gz' and arguments.vcf[-4:] != '.vcf':
            sys.exit('\nArgument --vcf accepts files ending in ".vcf.gz" or ".vcf".\n')

        elif arguments.genos or arguments.markers:
            sys.exit('\nExiting process...\nYou cannot use a vcf file together with a genotyping or markers file.\n')

        roh_detection(arguments, variable)

    if not arguments.vcf:

        if (arguments.genos and not arguments.markers) or (arguments.markers and
                                                           not arguments.genos):
            sys.exit('\nExiting process...\nArguments --markers and --genos should be used together\n')

        f2data_analysis(arguments, variable)
        plot_pairs(arguments)


if __name__ == '__main__':
    main()
