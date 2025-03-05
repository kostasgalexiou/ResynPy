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


def filter_hetero(content, hetero_threshold):
    # keep only individuals that have a heterozygosity level lower than the assigned argument
    new_content = list()
    disc_list = list()
    d_individuals = 0
    for g in content:
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
    pairs_genos_count = collections.Counter(sorted([''.join(sorted(p1 + p2)) for p1, p2 in zip(genos[i1_index][1:],
                                                                                               genos[i2_index][1:])]))

    hh_count = pairs_genos_count['HH']
    return sum([scores_dictionary[k] * v for k, v in pairs_genos_count.items()]), hh_count


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


def run_analysis(index_chunk, new_content, marker_file, ideal_score, scoresd,
                 invar, total_score_threshold, phased):
    disc_pair_list = list()
    discarded_pairs_invar = 0
    discarded_pairs_score = 0
    accepted_pairs = 0
    out_list = list()

    for i, (indv1, indv2) in enumerate(index_chunk):
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

        similarity_higher = '%.2f' % float((ideal_score - invar_counts) * 100 / ideal_score)
        similarity_lower = '%.2f' % float((ideal_score - (invar_counts + hhcount + het_only1 +
                                                          het_only2)) * 100 / ideal_score)
        accepted_pairs += 1

        if phased == 'yes':
            rec_count1 = recombination_count(m_file=marker_file, geno_string=new_content[indv1][1:])
            rec_count2 = recombination_count(m_file=marker_file, geno_string=new_content[indv2][1:])
            out_list.append(['|'.join((new_content[indv1][0], new_content[indv2][0])),
                             str(score_pair),
                             str(invar_counts),
                             heteros,
                             '-'.join((similarity_lower, similarity_higher)),
                             str(rec_count1), str(rec_count2), str(rec_count1 + rec_count2)])
        else:
            out_list.append(['|'.join((new_content[indv1][0], new_content[indv2][0])),
                             str(score_pair),
                             str(invar_counts), heteros,
                             '-'.join((similarity_lower, similarity_higher))])

    return out_list, accepted_pairs, disc_pair_list, discarded_pairs_score, discarded_pairs_invar
