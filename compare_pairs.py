# -*- coding: utf-8 -*-
"""
Created on {8/21/18}

@author: {kalexiou}
"""

import sys
import collections, itertools, warnings, os
from collections import OrderedDict as od

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

"""It takes ~6min to run the analysis for 2278 individuals, when a score filter of >0.5*ideal_score is applied 
(~2,590,000 comparisons)."""


def logfile(results_dir, string):
    with open('%s/PairComparison.log' % results_dir, 'a') as log:
        log.writelines(string)

    return


def scoresdict(score_file):
    scores_dict = od()

    with open(score_file) as f:
        for line in f:
            geno, score = line.rstrip('\n\r').split('\t')
            scores_dict[geno] = float(score)

    return scores_dict


def compare(p1, p2, dictionary):
    pair = '/'.join([p1, p2])
    reverse_pair = '/'.join([p2, p1])

    score = int()
    if pair in dictionary.keys():
        score = dictionary[pair]
    elif reverse_pair in dictionary.keys():
        score = dictionary[reverse_pair]

    return score


def filter_hetero(genotypes, hetero_threshold):
    # keep only individuals that have a heterozygosity level lower than the assigned argument
    new_content = list()
    disc_list = list()
    d_individuals = 0
    for g in genotypes:
        geno_with_info = [x for x in g if x != '-']
        hetero_count = geno_with_info.count('H')
        hetero_ratio = float(hetero_count) / float(len(geno_with_info))

        if hetero_ratio >= hetero_threshold:
            disc_list.append([g[0], 'heterozygosity: %.2f' % hetero_ratio])
            d_individuals += 1
            continue

        new_content.append(g)

    return new_content, disc_list, d_individuals


def hetero_average_and_invariable(genos, i1_index, i2_index, scores):
    hetero1_count = len([x for x in genos[i1_index][1:] if x == 'H'])
    wInfo1 = len([x for x in genos[i1_index][1:] if x != '-'])
    ratio1 = float(hetero1_count) / float(wInfo1)
    hetero2_count = len([x for x in genos[i2_index][1:] if x == 'H'])
    wInfo2 = len([x for x in genos[i2_index][1:] if x != '-'])
    ratio2 = float(hetero2_count) / float(wInfo2)
    average = float(sum([ratio1, ratio2])) / 2

    invar_genotypes = 0
    score_pair = 0
    invar_genotypes_list = list()
    for p1, p2 in zip(genos[i1_index][1:], genos[i2_index][1:]):
        if p1 == p2 and p1 != '-' and p1 != 'H':
            invar_genotypes += 1
            invar_genotypes_list.append(p1 + p2)
        score_pair += compare(p1=p1, p2=p2, dictionary=scoresdict(scores))

    return average, score_pair, invar_genotypes, invar_genotypes_list


def recombination_count(m_file, geno_string):
    chrom_list = list()
    with open(m_file) as f:
        for elem in f:
            chrom_list.append(elem.rstrip('\n\r').split('\t')[0])

    # get number of markers per chromosome
    count_dict = collections.OrderedDict(collections.Counter(sorted(chrom_list)))

    # generate strings of genotypes, grouping by number of markers per chromosome.
    chr_list = []
    marker_no = 0
    for i, key in enumerate(sorted(count_dict.keys())):

        if i == 0:
            chr_list.append(''.join(geno_string[:count_dict[key]]))
        else:
            start = marker_no
            end = marker_no + count_dict[key]
            chr_list.append(''.join(geno_string[start:end]))

        marker_no += count_dict[key]

    recomb_count = 0
    for item in chr_list:
        recomb_count += item.count('AHB')
        recomb_count += item.count('AB')
        recomb_count += item.count('BA')
        recomb_count += item.count('BHA')

    return recomb_count


def save_discarded(res_dir, prefix, list_of_discarded):
    discfile = '/'.join((res_dir, prefix + '_discarded_individuals-pairs.tab'))
    with open(discfile, 'w') as discarded:
        for elem in list_of_discarded:
            discarded.writelines('\t'.join(elem) + '\n')

    return


def main():
    infile = sys.argv[1]
    marker_file = sys.argv[2]
    results_folder = sys.argv[3]
    scoresfile = sys.argv[4]
    heterozygosity = float(sys.argv[5])
    invar = float(sys.argv[6])
    pref = sys.argv[7]
    phased = sys.argv[8]
    outfile = '/'.join((results_folder, pref + '_selected-pairs.tab'))

    with open(outfile, 'w') as out:

        if phased == 'yes':
            out.writelines('\t'.join(('#indiv_pair', 'score', 'invariable_pairs', 'average_heterozygosity',
                                      'recomb_no1', 'recomb_no2', 'sum_recomb\n',)))
        else:
            out.writelines('\t'.join(('#indiv_pair', 'score', 'invariable_pairs', 'average_heterozygosity\n')))

        # discarded_indvs = 0
        with open(infile, "rt") as f:

            f1 = f.readlines()
            content = [x.rstrip('\n\r').split('\t') for x in f1]

            # keep only individuals that have a heterozygosity level lower than the assigned argument
            new_content, disc_indv_list, d_individuals = filter_hetero(content, heterozygosity)

            logfile(results_folder, '\tComparison Statistics' + '\n')
            logfile(results_folder, '\t---------------------\n\n')
            logfile(results_folder, '\t\ttotal number of individuals  : ' + str(len(content[1:])) + '\n')
            logfile(results_folder, '\t\tindividuals not meeting\n'
                                    '\t\t  heterozygosity threshold   : ' +
                    str(len(content[1:]) - len(new_content[1:])) + '\n')
            logfile(results_folder, '\t\tindividuals to compare       : ' + str(len(new_content[1:])) + '\n')
            logfile(results_folder, '\t\tpairs to compare             : ' +
                    str((len(new_content[1:]) * (len(new_content[1:]) - 1)) / 2).split('.')[0] + '\n')

            index_pairs = list(itertools.combinations(range(1, len(new_content)), 2))

            with open(marker_file) as tmp:
                tmp1 = tmp.readlines()
                ideal_score: int = len(tmp1)

            disc_pair_list = list()
            discarded_pairs_invar = 0
            discarded_pairs_score = 0
            accepted_pairs = 0
            for indv1, indv2 in index_pairs:

                hetero_avg, score_pair, same_genotypes, same_genotypes_list = hetero_average_and_invariable(
                    genos=new_content, i1_index=indv1, i2_index=indv2, scores=scoresfile)

                # remove pairs with a score less that 80% of the ideal score or pairs whose ratio of
                # invariable sites is equal or higher than the filter_invar.
                if float(same_genotypes) / float(ideal_score) >= invar:
                    disc_pair_list.append([new_content[indv1][0] + '|' + new_content[indv2][0],
                                           'invariable genotype ratio: %.2f' %
                                           (float(same_genotypes) / float(ideal_score))])
                    discarded_pairs_invar += 1
                    continue

                elif score_pair < 0.8 * float(ideal_score):
                    disc_pair_list.append([new_content[indv1][0] + '|' + new_content[indv2][0],
                                           'pair score: %.2f' % score_pair])
                    discarded_pairs_score += 1
                    continue

                aa_counts = same_genotypes_list.count('AA')
                bb_counts = same_genotypes_list.count('BB')
                invar_counts = aa_counts + bb_counts

                accepted_pairs += 1
                if phased == 'yes':
                    rec_count1 = recombination_count(m_file=marker_file, geno_string=new_content[indv1][1:])
                    rec_count2 = recombination_count(m_file=marker_file, geno_string=new_content[indv2][1:])
                    out.writelines('\t'.join(('|'.join((new_content[indv1][0], new_content[indv2][0])),
                                              str(score_pair), str(invar_counts), '%.2f' % hetero_avg,
                                              str(rec_count1), str(rec_count2), str(rec_count1 + rec_count2) + '\n'
                                              )))
                else:
                    out.writelines('\t'.join(('|'.join((new_content[indv1][0], new_content[indv2][0])),
                                              str(score_pair), str(invar_counts), '%.2f' % hetero_avg + '\n')))

            logfile(results_folder, '\t\tdiscarded pairs due to score : ' +
                    str(discarded_pairs_score) + '\n')
            logfile(results_folder, '\t\tpairs not meeting invariable\n'
                                    '\t\t  site threshold             : ' +
                    str(discarded_pairs_invar) + '\n')
            logfile(results_folder, '\t\tACCEPTED PAIRS               : ' +
                    str(accepted_pairs) + '\n\n')

    # save discarded individuals and pairs into a separate file
    discarded_elements = disc_indv_list + disc_pair_list
    save_discarded(res_dir=results_folder, prefix=pref, list_of_discarded=discarded_elements)

    if sys.platform == "linux" or sys.platform == "linux2":
        command0 = 'export LC_NUMERIC=en_US.utf-8'
        os.system(command0)

        if phased == 'yes':
            command1a = ' '.join(('head -1', outfile, '> a'))
            command1b = ' '.join(('tail -n+2', outfile, '| sort -g -k2,2nr -k7,7nr >', 'b'))
            command1c = 'cat a b > x'
        elif phased == 'no':
            command1a = ' '.join(('head -1', outfile, '> a'))
            command1b = ' '.join(('tail -n+2', outfile, '| sort -g -k2,2nr >', 'b'))
            command1c = 'cat a b > x'

        os.system(command1a)
        os.system(command1b)
        os.system(command1c)

        command2 = ' '.join(('mv x', outfile))
        os.system(command2)
        os.system('rm a b')

    elif sys.platform == 'win32':
        print('\nYou are using a win32 operating system!\n\n')
        print('\tPlease remember to order manually file %s+"_selected-pairs.tab", in order to have the best '
              'pairs on the top of the list.\n' % pref)


if __name__ == '__main__':
    main()
