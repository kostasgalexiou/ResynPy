# -*- coding: utf-8 -*-
"""
Created on {8/21/18}

@author: {kalexiou}
"""

import sys, collections, itertools, warnings, os
import os.path as op
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


def recombination_count(m_file, geno_string):
    chrom_list = list()
    with open(m_file) as f:
        for elem in f:
            chrom_list.append(elem.rstrip('\n\r').split('\t')[1])

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


def main():
    infile = sys.argv[1]
    marker_file = sys.argv[2]
    results_folder = sys.argv[3]
    scoresfile = sys.argv[4]
    heterozygosity = float(sys.argv[5])
    filter_invar = float(sys.argv[6])
    phased = sys.argv[7]

    outfile = '/'.join((results_folder, 'individual_genetic_distance.tab'))

    with open(outfile, 'w') as out:

        with open(results_folder+'/discarded_individuals_and_pairs.tab', 'w') as discarded:

            if phased == 'yes':
                out.writelines('\t'.join(('#indiv_pair', 'score', 'AA pairs', 'BB pairs', 'recomb_no1', 'recomb_no2',
                                          'sum_recomb\n',)))
            else:
                out.writelines('\t'.join(('#indiv_pair', 'score', 'AA pairs', 'BB pairs\n')))

            with open(infile, "rt") as f:

                f1 = f.readlines()
                content = [x.rstrip('\n\r').split('\t') for x in f1]

                # keep only individuals that have a heterozygosity level lower than the assigned argument
                new_content = list()
                for elem in content:
                    new_elem = ''
                    geno_with_info = [x for x in elem if x != '-']
                    hetero_count = geno_with_info.count('H')
                    hetero_ratio = float(hetero_count) / float(len(geno_with_info))

                    if hetero_ratio >= heterozygosity:
                        discarded.writelines('\t'.join((elem[0], 'heterozygosity: %.2f' % hetero_ratio+'\n')))
                        continue

                    new_elem += '|'.join(elem)
                    new_content.append(new_elem.split('|'))

                index_pairs = list(itertools.combinations(range(1, len(new_content)), 2))

                with open(marker_file) as tmp:
                    tmp1 = tmp.readlines()
                    ideal_score: int = len(tmp1)

                for indv1, indv2 in index_pairs:

                    same_genotypes = 0
                    score_pair = 0
                    same_genotypes_list = list()
                    for p1, p2 in zip(new_content[indv1][1:], new_content[indv2][1:]):
                        if p1 == p2 and p1 != '-' and p1 != 'H':
                            same_genotypes += 1
                            same_genotypes_list.append(p1+p2)
                        score_pair += compare(p1=p1, p2=p2, dictionary=scoresdict(scoresfile))

                    # remove pairs with a score less that 80% of the ideal score or pairs whose percentage of
                    # identical pairs are equal or higher than the filter_invar.
                    if float(same_genotypes) / float(ideal_score) >= filter_invar:
                        discarded.writelines('\t'.join((new_content[indv1][0]+'|'+new_content[indv2][0],
                                                        'invariable genotype ratio: %.2f' % (float(same_genotypes) /
                                                                                             float(ideal_score))+'\n')))
                        continue

                    elif score_pair < 0.8 * float(ideal_score):
                        discarded.writelines('\t'.join((new_content[indv1][0]+'|'+new_content[indv2][0],
                                                        'pair score: %.2f' % score_pair+'\n')))
                        continue

                    aa_counts = same_genotypes_list.count('AA')
                    bb_counts = same_genotypes_list.count('BB')

                    if phased == 'yes':
                        rec_count1 = recombination_count(marker_file, new_content[indv1][1:])
                        rec_count2 = recombination_count(marker_file, new_content[indv2][1:])
                        out.writelines('\t'.join(('|'.join((new_content[indv1][0], new_content[indv2][0])),
                                                  str(score_pair), str(aa_counts), str(bb_counts),
                                                  str(rec_count1), str(rec_count2), str(rec_count1 + rec_count2) + '\n'
                                                  )))
                    else:
                        out.writelines('\t'.join(('|'.join((new_content[indv1][0], new_content[indv2][0])),
                                                  str(score_pair), str(aa_counts), str(bb_counts) + '\n')))

    if sys.platform == "linux" or sys.platform == "linux2":
        if phased == 'yes':
            command1 = ' '.join(('sort -g -k2,2r -k7,7nr', outfile, '>', 'x'))
        elif phased == 'no':
            command1 = ' '.join(('sort -g -k2,2r', outfile, '>', 'x'))

        os.system(command1)
        logfile(results_folder, '\n' + '\t' + command1 + '\n')
        command2 = ' '.join(('mv x', outfile))
        os.system(command2)
        logfile(results_folder, '\t' + command2 + '\n\n')

    elif sys.platform == 'win32':
        print('\n\tPlease order manually file "individual_genetic_distance.tab" based on the score column\n')





if __name__ == '__main__':
    main()
