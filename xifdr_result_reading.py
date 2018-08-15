import pandas as pd


def build_dfs_of_int_betw(f):
    """
    :returns df_between and df_internal, only containing TT links
    """
    df = pd.read_csv(f)
    df = df[df['isTT']]
    df_internal = df[df['fdrGroup'].str.contains('Within|[Ii]nternal')]
    df_between = df[df['fdrGroup'].str.contains('[Bb]etween')]
    if not len(df) == len(df_internal) + len(df_between):
        raise ValueError("Input Dataframe must have same number of 'TT' rows as combined output Dataframes")
    return df_between, df_internal


def scores_from_filedict(dct_f):
    def f(key):
        _, n1, _, n2 = key.split("_")
        return float(n1) + float(n2 ) /10000
    df_scores = pd.DataFrame()
    for db in sorted(dct_f, key=f):
        files = dct_f[db]
        #         print db
        for i, f in enumerate(files):
            idx = db
            replicate = i+ 1
            var_int = "run {}: internal TT".format(replicate)
            var_bet = "run {}: between TT".format(replicate)

            df_bet, df_int = build_dfs_of_int_betw(f)

            df_sc_bet = pd.DataFrame(df_bet["Score"])
            df_sc_bet["id"] = idx
            df_sc_bet["variable"] = var_bet

            df_sc_int = pd.DataFrame(df_int["Score"])
            df_sc_int["id"] = idx
            df_sc_int["variable"] = var_int
            df_scores = df_scores.append([df_sc_int, df_sc_bet], ignore_index=True)
    df_scores["link type"] = df_scores["variable"].map(lambda x: x.split(": ")[1])
    return df_scores


def minimums_from_scores_df(df_scores):
    df = pd.DataFrame()
    for db_comp in df_scores["id"].unique():
        df_runs = df_scores[df_scores["id"].str.match(db_comp)]
        for run in df_runs["variable"].unique():
            df_r = df_runs[df_runs["variable"].str.match(run)]
            df.loc[db_comp, run] = df_r["Score"].min()
    return df


def build_long_min_df(df_scores):
    df_wide = minimums_from_scores_df(df_scores)
    df_wide["id"] = df_wide.index
    df_sc_min_long = df_wide.melt(id_vars="id", value_name="Score").dropna()

    def f(string):
        return string.split(" ")[2]

    df_sc_min_long["link type"] = df_sc_min_long["variable"].map(f)

    def f(key):
        _, n1, _, n2 = key.split("_")
        return n1 + "\n" + n2

    df_sc_min_long["DB composition"] = df_sc_min_long["id"].map(f)

    def f(key):
        _, n1, _, n2 = key.split("_")
        return int(n1) + int(n2)

    df_sc_min_long["total DB size"] = df_sc_min_long["id"].map(f)
    return df_sc_min_long
