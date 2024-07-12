# -*- coding: utf-8 -*-
"""
Created on {8/21/18}

@author: {kalexiou}
"""

import collections
import itertools
import time
import warnings

import pandas as pd

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")


def logfile(results_dir, string):
    with open('%s/PairComparison.log' % results_dir, 'a') as log:
        log.writelines(string)

    return


def filter_hetero(genotypes, hetero_threshold):
    # keep only individuals that have a heterozygosity level lower than the assigned argument
    new_content = list()
    disc_list = list()
    d_individuals = 0
    for g in genotypes:
        geno_with_info = [x for x in g if x != '-']
        hetero_count = geno_with_info.count('H')
        hetero_ratio = float(hetero_count) / float(len(geno_with_info))

        if hetero_ratio >= float(hetero_threshold):
            disc_list.append([g[0], 'heterozygosity: %.2f' % hetero_ratio])
            d_individuals += 1
            continue

        new_content.append(g)

    return new_content, disc_list, d_individuals


def hetero_average_and_invariable(genos, i1_index, i2_index):

    hetero1_count = genos[i1_index][1:].count('H')
    hetero2_count = genos[i2_index][1:].count('H')
    both_heteros_count = len([i for i, (a, b) in enumerate(zip(genos[i1_index][1:], genos[i2_index][1:])) if a == b
                              and a == 'H'])
    only1_hetero = hetero1_count - both_heteros_count
    only2_hetero = hetero2_count - both_heteros_count

    winfo1 = len([x for x in genos[i1_index][1:] if x != '-'])
    ratio1 = float(hetero1_count) / float(winfo1)

    winfo2 = len([x for x in genos[i2_index][1:] if x != '-'])
    ratio2 = float(hetero2_count) / float(winfo2)

    heterozygosities = '|'.join(('%.2f' % ratio1, '%.2f' % ratio2))

    invar_genotypes_list = []
    for p1, p2 in zip(genos[i1_index][1:], genos[i2_index][1:]):
        if p1 == p2 and p1 != '-' and p1 != 'H':
            invar_genotypes_list.append(p1 + p2)

    return heterozygosities, len(invar_genotypes_list), invar_genotypes_list, only1_hetero, only2_hetero


def total_score(genos, i1_index, i2_index, scores_dictionary):

    pairs_genos_count = collections.Counter(sorted([''.join(sorted(p1+p2)) for p1, p2 in zip(genos[i1_index][1:],
                                                                                             genos[i2_index][1:])]))

    hh_count = pairs_genos_count['HH']
    return sum([scores_dictionary[k]*v for k, v in pairs_genos_count.items()]), hh_count


def recombination_count(m_file, geno_string):
    chrom_list = list()
    with open(m_file) as f:
        for elem in f:
            chrom_list.append(elem.rstrip('\n\r').split('\t')[0])

    # get the number of markers per chromosome
    count_dict = collections.OrderedDict(collections.Counter(sorted(chrom_list)))

    # generate strings of genotypes, grouping by the number of markers per chromosome.
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


def run_analysis(infile, marker_file, results_folder, scoresd, heterozygosity, invar, total_score_threshold, pref,
                 phased):

    outfile = '/'.join((results_folder, pref + '_selected-pairs.tab'))

    with open(outfile, 'w') as out:

        if phased == 'yes':
            out.writelines('\t'.join(('#indiv_pair', 'score', 'invariable_pairs', 'hetero1|hetero2',
                                      'recomb_no1', 'recomb_no2', 'sum_recomb', 'similarity_to_hybrid\n')))
        else:
            out.writelines('\t'.join(('#indiv_pair', 'score', 'invariable_pairs', 'hetero1|hetero2',
                                      'similarity_to_hybrid\n')))

        with open(infile, "rt") as f:

            f1 = f.readlines()
            content = [x.rstrip('\n\r').split('\t') for x in f1]

            # keep only individuals that have a heterozygosity level lower than the assigned argument
            new_content, disc_indv_list, d_individuals = filter_hetero(content, heterozygosity)

            logfile(results_folder, '\tComparison Statistics' + '\n')
            logfile(results_folder, '\t---------------------\n\n')

            with open(marker_file) as tmp:
                tmp1 = tmp.readlines()
                ideal_score: int = len(tmp1)

            logfile(results_folder, '\t\ttotal number of markers      : ' + str(ideal_score) + '\n')
            logfile(results_folder, '\t\ttotal number of individuals  : ' + str(len(content[1:])) + '\n')
            logfile(results_folder, '\t\tindividuals not meeting\n'
                                    '\t\t  heterozygosity threshold   : ' +
                    str(len(content[1:]) - len(new_content[1:])) + '\n')
            logfile(results_folder, '\t\tindividuals to compare       : ' + str(len(new_content[1:])) + '\n')
            logfile(results_folder, '\t\tpairs to compare             : ' +
                    str((len(new_content[1:]) * (len(new_content[1:]) - 1)) / 2).split('.')[0] + '\n')

            index_pairs = list(itertools.combinations(range(1, len(new_content)), 2))

            disc_pair_list = list()
            discarded_pairs_invar = 0
            discarded_pairs_score = 0
            accepted_pairs = 0
            start = time.perf_counter()

            for i, (indv1, indv2) in enumerate(index_pairs):
                if i % 10000 == 0:
                    if i == 0:
                        print('I have started the comparison of {0} individual pairs. '
                              'Current time --> {1}'.format(len(index_pairs), time.strftime('%X')))
                    else:
                        stop = time.perf_counter()
                        difference: str = '%.2f' % (float(stop-start))
                        print('I have analyzed {0} out of {1} pairs in {2} seconds'.format(i, len(index_pairs),
                                                                                           difference))
                        start = time.perf_counter()

                heteros, same_genotypes, same_genotypes_list, het_only1, het_only2 = hetero_average_and_invariable(
                    genos=new_content, i1_index=indv1, i2_index=indv2)

                # filter for invariable sites and total score
                if float(same_genotypes) / float(ideal_score) >= float(invar):
                    disc_pair_list.append([new_content[indv1][0] + '|' + new_content[indv2][0],
                                           'invariable genotype ratio: %.2f' %
                                           (float(same_genotypes) / float(ideal_score))])
                    discarded_pairs_invar += 1
                    continue

                else:

                    score_pair, hhcount = total_score(genos=new_content, i1_index=indv1, i2_index=indv2,
                                             scores_dictionary=scoresd)
                    if score_pair < float(total_score_threshold) * float(ideal_score):
                        disc_pair_list.append([new_content[indv1][0] + '|' + new_content[indv2][0],
                                               'pair score: %.2f' % score_pair])
                        discarded_pairs_score += 1
                        continue

                aa_counts = same_genotypes_list.count('AA')
                bb_counts = same_genotypes_list.count('BB')
                invar_counts = aa_counts + bb_counts

                similarity_higher = '%.2f' % float((ideal_score - invar_counts)*100/ideal_score)
                similarity_lower = '%.2f' % float((ideal_score - (invar_counts + hhcount + het_only1 +
                                                                  het_only2))*100/ideal_score)

                accepted_pairs += 1
                if phased == 'yes':
                    rec_count1 = recombination_count(m_file=marker_file, geno_string=new_content[indv1][1:])
                    rec_count2 = recombination_count(m_file=marker_file, geno_string=new_content[indv2][1:])
                    out.writelines('\t'.join(('|'.join((new_content[indv1][0], new_content[indv2][0])),
                                              str(score_pair), str(invar_counts), heteros,
                                              str(rec_count1), str(rec_count2), str(rec_count1 + rec_count2),
                                              '-'.join((similarity_lower, similarity_higher)) + '\n'
                                              )))
                else:
                    out.writelines('\t'.join(('|'.join((new_content[indv1][0], new_content[indv2][0])),
                                              str(score_pair), str(invar_counts), heteros,
                                              '-'.join((similarity_lower, similarity_higher)) + '\n')))

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

    # order file with selected pairs
    a = pd.read_table(outfile, header=0, sep='\t')
    b = a.sort_values(by=["score", "sum_recomb"], ascending=[False, True])
    b.to_csv(outfile, index=False, sep='\t')

    return
