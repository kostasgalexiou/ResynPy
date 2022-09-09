#!/usr/bin/python3

"""
@filename: plot_pair_genotypes.py
@author: kalexiou
@date: 5/9/22
"""

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import sys
from collections import OrderedDict as od
from collections import Counter
from matplotlib.lines import Line2D

pairs = sys.argv[1]  # compare_pairs.py output
genotypes = sys.argv[2]  # genotyping data
markers = sys.argv[3]
figname_prefix = sys.argv[4]


def get_total_markers():
    # get number of markers and number of markers per chromosome
    with open(markers) as m:
        m1 = m.readlines()
        total_markers = len(m1)

        chr_list = list()
        for elem in m1:
            chr_list.append(elem.split('\t')[0])

        chr_list_unique = sorted(list(set(chr_list)))
        number_per_chr = od(sorted(Counter(chr_list).items()))

        g_indices = list()
        gap = 0
        for m in list(number_per_chr.values())[:-1]:
            gap += m
            g_indices.append(gap)

        gap_indices = [x - 1 for x in g_indices]  # 0-based

    return total_markers, gap_indices, chr_list_unique


def indv2genos():  # get the correspondence of individual to its genotype
    dictionary = od()
    with open(genotypes) as f:
        f_list = f.readlines()
        for i in f_list[1:]:
            individual = i.rstrip('\n\r').split('\t')[0]
            genos = ','.join(i.rstrip('\n\r').split('\t')[1:])
            dictionary[individual] = genos

    return dictionary


def get_indv_colors(ind_genos, white_pos):
    ind_colors = list()
    color = ''
    for i, geno in enumerate(ind_genos):

        if i in white_pos:
            if geno == 'A':
                ind_colors.append('salmon')
            elif geno == 'B':
                ind_colors.append('palegreen')
            elif geno == 'H':
                ind_colors.append('yellow')
            elif geno == '-':
                ind_colors.append('lightgrey')

            ind_colors.append('white')

        else:
            if geno == 'A':
                color = 'salmon'
            elif geno == 'B':
                color = 'palegreen'
            elif geno == 'H':
                color = 'yellow'
            elif geno == '-':
                color = 'lightgrey'
            ind_colors.append(color)

    return ind_colors


def generate_graph(indv_dict):

    all_markers, white, chromosomes = get_total_markers()

    with open(pairs) as f:

        f1 = f.readlines()

        genotype_color_lists = list()
        line_names = list()
        for i, n in enumerate(f1[1:11]):

            pair = n.rstrip('\n\r').split('\t')[0].split('|')
            line_names.append(pair)

            ind1_genos = indv_dict[pair[0]].split(',')
            ind2_genos = indv_dict[pair[1]].split(',')

            ind1_colors = get_indv_colors(ind1_genos, white)
            ind2_colors = get_indv_colors(ind2_genos, white)
            genotype_color_lists.append([ind1_colors, ind2_colors])

    fig, ax = plt.subplots()

    ax.set_ylim(0, all_markers + int(all_markers/4))  # y-axis limits
    ax.set_xticks(range(0, all_markers + len(white) + 5, 5))  # x-axis limits
    plt.rc('font', size=4)

    range_of_ys = range(int((all_markers + 10) / 10), all_markers + 10, int((all_markers + 10) / 10))[:10]  # start,
    # stop, step (7, 72, 7); first 10 values only

    for y, color_list_indvs, l_name in zip(range_of_ys, genotype_color_lists, line_names):

        if len(ax.set_xticks(range(0, all_markers + len(white) + 5, 5))) < 15:
            ax.text(-5, y-0.5, l_name[0])
            ax.text(-5, y+max(plt.ylim(0, all_markers + int(all_markers/4))) / 60, l_name[1])
        else:
            ax.text(-10, y - 0.5, l_name[0])
            ax.text(-10, y + max(plt.ylim(0, all_markers + int(all_markers / 4))) / 60, l_name[1])

        r = -1
        c_index = 0
        for c1, c2, x in zip(color_list_indvs[0], color_list_indvs[1], range(1, all_markers + len(white)+1, 1)):
            r += 1
            if x == 1 and y == max(range_of_ys):  # -int((all_markers + 20) / 10):
                ax.text(x-0.7, y+(max(plt.ylim(0, all_markers + int(all_markers/4))) / 60)*2.5, chromosomes[0])

            if c1 == 'white' and y == max(range_of_ys):  # -int((all_markers + 20) / 10):
                c_index += 1
                ax.text(range(1, all_markers + len(white)+1, 1)[r+1] - 0.7, y + (max(plt.ylim(0, all_markers + int(all_markers/4))) / 60)*2.5, chromosomes[c_index])
                continue

            if c1 == 'white':
                continue

            ax.add_patch(Rectangle((x-0.7, y-0.5), 1, max(plt.ylim(0, all_markers + int(all_markers/4))) / 60,
                                   ec='darkgrey', linewidth=0.2, color=c1))
            ax.add_patch(Rectangle((x-0.7, y + max(plt.ylim(0, all_markers + int(all_markers/4))) / 60), 1,
                                   max(plt.ylim(0, all_markers + int(all_markers/4))) / 60, ec='darkgrey', linewidth=0.2,
                                   color=c2))

            if y < range_of_ys[-1]:
                ax.axhline(y+int(int((all_markers + 10) / 10)/2)+1, xmin=0, xmax=1, lw=0.5, color='darkgrey')

    custom_lines = [Line2D([0], [0], color='salmon', lw=4),
                    Line2D([0], [0], color='palegreen', lw=4),
                    Line2D([0], [0], color='yellow', lw=4),
                    Line2D([0], [0], color='lightgrey', lw=4)]

    ax.legend(custom_lines, ['A', 'B', 'H', 'missing'], loc='lower center', ncol=4, frameon=False)

    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    for spine in ['top', 'right', 'left', 'bottom']:
        ax.spines[spine].set_visible(False)

    ax.set_facecolor('azure')
    ax.set_title('Top 10 of selected pairs of invididuals', size=6)

    my_dpi = 96
    plt.savefig('%s.png' % figname_prefix, dpi=my_dpi * 20)
    plt.savefig('%s.pdf' % figname_prefix, dpi=my_dpi * 20)
    # plt.show()


if __name__ == '__main__':

    ind_dict = indv2genos()
    generate_graph(indv_dict=ind_dict)
