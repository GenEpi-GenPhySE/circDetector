#!/usr/bin/env python3

# Imports
from collections import defaultdict
import operator, pandas, re
import pybedtools as pybed
import natsort
from tqdm import tqdm
import pandas as pd


def intersection(lst1, lst2): 
    return list(set(lst1) & set(lst2)) 


names=["cow", "pig"]
df_ortholoy = pd.read_table("1796_orthoBovPig-exonicCircRNA.txt", sep="\t", index_col=0, names=names).fillna(" ")
df_all_circ = pd.read_table("cirRNAcounts.tsv", sep="\t", index_col=0).fillna(" ")

circ_ortho_pig = df_ortholoy["pig"].index.to_list()
circ_ortho_cow = df_ortholoy["cow"].index.to_list()

circ_ortho_all = []
circ_ortho_all.extend(circ_ortho_pig)
circ_ortho_all.extend(circ_ortho_cow)

# circ_not_found = []
# for circ in circ_ortho_all:
#     if circ not in df_all_circ.index:
#         circ_not_found.append(circ)
# print(len(circ_not_found))

df_final = df_all_circ[df_all_circ.index.isin(circ_ortho_all)]
df_final.to_csv("cirRNAcounts_ortholog.tsv", sep="\t", index=True)