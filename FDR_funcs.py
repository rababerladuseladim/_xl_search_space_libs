import pandas as pd

g_enzyme_regex_dict = {
    "trypsin": r"[KR][^P]+"
}

def split_mod_and_unmod_peptides(df_pep, modifications=["Mox", "bs3"]):
    regex_pattern = "|".join(modifications)
    s_bool = (df_pep["Peptide1"].str.contains(regex_pattern) | df_pep["Peptide2"].str.contains(regex_pattern))
    df_mod = df_pep[s_bool]
    df_umod = df_pep[~s_bool]
    return df_umod, df_mod


def calc_FDR(df):
    """(TD - DD)/TT"""
    fdr = (df["isTD"].sum() - df["isDD"].sum()) / float(df["isTT"].sum())
    return fdr


def FDR_for_mod_subsets(df, **kwargs):
    """
    calculates FDR for modified/unmodified jpeptides separately
    kwargs are handed over to split_mod_and_unmod_peptides
    """
    df_super = df
    df_umod, df_mod = split_mod_and_unmod_peptides(df, **kwargs)
    col_no = "no. [TT]"
    col_FDR = "FDR [%]"
    df_ret = pd.DataFrame(columns=[col_FDR, col_no])
    for df, db in zip([df_super, df_umod, df_mod], ["total", "unmodified", "modified"]):
        df_ret.loc[db, col_FDR] = calc_FDR(df)*100
        df_ret.loc[db, col_no] = df["isTT"].sum()
    return df_ret


def count_missed_cleavage(ser_of_peptide_str, regex_cleav):
    """
    return series with number of regex occurence in each peptide string. Before comparison, peptide strings are stripped
    from any non-amino-acid charecters and strings before and after cleavage sites are removed.
    For empty strings, the return value is 0.
    :param ser_of_peptide_str:
    :param regex_cleav:
    :return:
    """
    # remove everything before and after cleavage sites
    s = ser_of_peptide_str.str.replace(r".*\.(.*)\..*", r"\1")
    # keep only amino acid letters
    s = s.str.replace(r"[^A-Z]", "")
    c = s.str.count(regex_cleav)
    # linear matches don't contain peptide2, the resulting nans are replaced to keep the series in int format
    c = c.fillna(0)
    return c.astype(int)
