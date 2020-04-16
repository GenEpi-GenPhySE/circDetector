#!/usr/bin/env python3
#usage: python3 scripts/plots_log_final_out.py -i mapping_stats_samples.tsv -o output.test

#imports:
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from itertools import groupby
import os
import argparse
import re
import seaborn as sns

def read_file(file):
    """ Read the sample file containing path files and return a pandas dataframe"""
    sample = pd.read_csv(file, sep='\t', dtype=str)
    return sample


def nb_read():
    data_table = pd.DataFrame({'Individual':['cow-testis-neonat1']*2 + ['cow-liver-B004']*2 + ['pig-testis-31']*2 + ['pig-liver-old2']*2,
                               'Read':(['R1'] + ['R2'])*4,
                               'Number of input reads':[73041056,73041056,73041056,73041056,73041056,73041056,73041056,73041056],
                               'Uniquely mapped reads number':[37438675,37438675,37438675,37438675,37438675,37438675,37438675,37438675],
                               'Number of reads mapped to multiple loci':[33825194,33825194,33825194,33825194,33825194,33825194,33825194,33825194],
                               'Number of reads mapped to too many loci':[12158,12158,12158,12158,12158,12158,12158,12158],
                               'Number of chimeric reads':[735392,735392,735392,735392,735392,735392,735392,735392]
                               })
    return data_table


def test_table2():
    data_table = pd.DataFrame({'Species':['Pig-liver']*2 + ['Pig-testis']*2 + ['Cow-testis']*2 + ['Cow-liver']*2,
                               'Read':(['R1'] + ['R2'])*4,
                               'Total number of splices':[6511483,6511483,6511483,6511483,6511483,6511483,6511483,6511483],
                               'Number of splices : Annotated':[5751924,5751924,5751924,5751924,5751924,5751924,5751924,5751924],
                               'Number of splices: Non-canonical':[42774,42774,42774,42774,42774,42774,42774,42774]
                               })
    return data_table


def test_table3():
    data_table = pd.DataFrame({'Species':['Pig-liver']*2 + ['Pig-testis']*2 + ['Cow-testis']*2 + ['Cow-liver']*2,
                               'Read':(['R1'] + ['R2'])*4,
                               'GT/AG':[6411237,6411237,6411237,6411237,6411237,6411237,6411237,6411237],
                               'GC/AG':[51218,51218,51218,51218,51218,51218,51218,51218],
                               'AT/AC':[6254,6254,6254,6254,6254,6254,6254,6254]
                               })
    return data_table


def test_table4():
    data_table = pd.DataFrame({'Species':['Pig-liver']*2 + ['Pig-testis']*2 + ['Cow-testis']*2 + ['Cow-liver']*2,
                               'Read':(['R1'] + ['R2'])*4,
                               'Size_read':[113,118,113,118,113,118,113,118]
                               })
    return data_table


def add_line(ax, xpos, ypos):
    line = plt.Line2D([xpos, xpos], [ypos + .1, ypos],
                      transform=ax.transAxes, color='black')
    line.set_clip_on(False)
    ax.add_line(line)

def label_len(my_index,level):
    labels = my_index.get_level_values(level)
    return [(k, sum(1 for i in g)) for k,g in groupby(labels)]

def label_group_bar_table(ax, df):
    ypos = -.1
    scale = 1./df.index.size
    for level in range(df.index.nlevels)[::-1]:
        pos = 0
        for label, rpos in label_len(df.index,level):
            lxpos = (pos + .5 * rpos)*scale
            ax.text(lxpos, ypos, label, ha='center', transform=ax.transAxes)
            add_line(ax, pos*scale, ypos)
            pos += rpos
        add_line(ax, pos*scale , ypos)
        ypos -= .1


# STATS:

def main():

    # df = pd.DataFrame(columns=["App","Feature1", "Feature2","Feature3",
    #                            "Feature4","Feature5",
    #                            "Feature6","Feature7","Feature8"], 
    #                   data=[["SHA",0,0,1,1,1,0,1,0],
    #                         ["LHA",1,0,1,1,0,1,1,0],
    #                         ["DRA",0,0,0,0,0,0,1,0],
    #                         ["FRA",1,0,1,1,1,0,1,1],
    #                         ["BRU",0,0,1,0,1,0,0,0],
    #                         ["PAR",0,1,1,1,1,0,1,0],
    #                         ["AER",0,0,1,1,0,1,1,0],
    #                         ["SHE",0,0,0,1,0,0,1,0]])

    # sns.set()
    # df.set_index('App').T.plot(kind='bar', stacked=True)

    # plt.figure(figsize=(12,8))
    # ax = sns.barplot(data=df, palette=sns.color_palette("GnBu", 10))
    # plt.xticks(rotation='vertical')
    # ax.grid(b=True, which='major', color='#d3d3d3', linewidth=1.0)
    # ax.grid(b=True, which='minor', color='#d3d3d3', linewidth=0.5)
    # plt.show()

    # df_init = read_file(args.input_file)
    # columns_name = list(df_init)
    # print(columns_name)

    # df = pd.DataFrame(columns=columns_name, 
    #                   data=[["SHA",0,0,1,1,1,0,1,0],
    #                         ["LHA",1,0,1,1,0,1,1,0],
    #                         ["DRA",0,0,0,0,0,0,1,0],
    #                         ["FRA",1,0,1,1,1,0,1,1],
    #                         ["BRU",0,0,1,0,1,0,0,0],
    #                         ["PAR",0,1,1,1,1,0,1,0],
    #                         ["AER",0,0,1,1,0,1,1,0],
    #                         ["SHE",0,0,0,1,0,0,1,0]])
    

    df = nb_read().groupby(['Individual','Read']).sum()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    df.plot(kind='bar',stacked=True,ax=fig.gca())
    #Below 3 lines remove default labels
    labels = ['' for item in ax.get_xticklabels()]
    ax.set_xticklabels(labels)
    ax.set_xlabel('')
    label_group_bar_table(ax, df)
    fig.subplots_adjust(bottom=.1*df.index.nlevels)
    plt.show()


def parse_arguments():
    parser = argparse.ArgumentParser(description='Sample file')
    parser.add_argument('-i', '--input_file',
                        required=True, help='Sample file')
    parser.add_argument('-o', '--output_file',
                        required=True, help='Sample file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_arguments()
    main()