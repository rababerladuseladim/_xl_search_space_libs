#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: henning
"""

import os
import pandas as pd

g_dct_occm_gene_to_uniprot = {
    'ORC1': 'P54784',
    'ORC2': 'P32833',
    'ORC3': 'P54790',
    'ORC4': 'P54791',
    'ORC5': 'P50874',
    'ORC6': 'P38826',
    'CDC6': 'P09119',
    'MCM2': 'P29469',
    'MCM3': 'P24279',
    'MCM4': 'P30665',
    'MCM5': 'P29496',
    'MCM6': 'P53091',
    'MCM7': 'P38132',
    'CDT1': 'P47112'
}


def occm_dict_from_fasta(fasta):
    dct_gene_uniprot_id = dict()
    with open(fasta) as f:
        for l in f.readlines():
            if l.startswith(">"):
                _, uniprot_id, fasta_desc = l.split("|")
                # gene = fasta_desc.split("GN=")[-1].split(" ")[0]
                gene = fasta_desc.split("_")[0]
                if gene not in dct_gene_uniprot_id:
                    dct_gene_uniprot_id[gene] = uniprot_id
                else:
                    raise AttributeError("'{}' appears more than once in input fasta!")
    return dct_gene_uniprot_id

# f = r"/home/henning/mnt/xitu/Data/Input/170322_OCCM_random_DB_test/fastas/occm_scerevisiae-180726/uniprot-yourlist%252525252525253AM201807268A530B6CA0138AFAA6D2B97CE8C2A92422BAE7U.fasta"


def rename_proteins(df_xifdr, dct_gene_uniprot=g_dct_occm_gene_to_uniprot):
    def rename_protein(protein_entry):
        p_ret = ""
        for i, p in enumerate(protein_entry.split(";")):
            if i > 0:
                p_ret += ";"
            p = p.upper()
            if "DECOY:" in p:
                p = p.split("DECOY:")[-1]
                p_ret += "DECOY:"
            if p in dct_gene_uniprot:
                p = dct_gene_uniprot[p]
            p_ret += p
        return p_ret

    for c in ["Protein1", "Protein2"]:
        df_xifdr[c] = df_xifdr[c].map(rename_protein)
    return df_xifdr


def do_not_execute_me():
    f = r"/home/henning/mnt/xitu/Data/Input/Modelling/OCCM/input_crosslinks/20180725-linkfdr_5/joined_links.csv"
    df = pd.read_csv(f)
    df_renamed = rename_proteins(df_xifdr=df)
    f, ext = os.path.splitext(f)
    f_out = f + "-renamed_prots" + ext
    df_renamed.to_csv(f_out)


