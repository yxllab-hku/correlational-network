import pandas as pd


def convert_float_to_int(_df: pd.DataFrame):
    float_col = _df.select_dtypes(include=['float'])
    float_col = float_col.astype('int')
    return pd.concat([_df.drop(float_col.columns, axis=1), float_col], axis=1)


def main():
    print("reading tables")
    genome = pd.read_excel("supplementary file 1.xlsx", sheet_name="Sheet10")["genome"]
    df_pre = pd.read_excel("supplementary file 1.xlsx", sheet_name="Sheet8")
    df_peptidase = pd.read_csv("23777967_protease_cluster.csv")

    # filter big clusters
    df_pre = df_pre[df_pre["number of members in cluster"] >= 10]
    df_pre['precursor group'] = 'pre_' + df_pre['precursor group'].apply(str)
    df_peptidase = df_peptidase[df_peptidase["number of members in cluster"] >= 100]
    df_peptidase['cluster No.'] = 'peptidase_' + df_peptidase['cluster No.'].apply(str)
    df_peptidase["genome"] = df_peptidase["mem"].str.replace(r"____.*$", "")

    # build count table
    df_occurrence = pd.DataFrame(genome)
    print("building count table pre")
    df_pre_piv = df_pre.pivot_table(values='sequence', index=['genome'], columns=['precursor group'], aggfunc='count')
    df_pre_piv.reset_index(inplace=True)
    df_occurrence = df_occurrence.merge(df_pre_piv, how='left')
    print("building count table pep")
    df_pep_piv = df_peptidase.pivot_table(values='mem', index=['genome'], columns=['cluster No.'], aggfunc='count')
    df_pep_piv.reset_index(inplace=True)
    df_occurrence = df_occurrence.merge(df_pep_piv, how='left')
    print("converting dtype")
    df_occurrence.fillna(0, inplace=True)
    df_occurrence = convert_float_to_int(df_occurrence)

    genome1 = pd.read_excel("supplementary file 1.xlsx", sheet_name="Sheet10").drop(['RefSeq accession', 'name'],
                                                                                    axis=1)
    out = pd.merge(genome1, df_occurrence)
    print("writing df_count.csv")
    out.to_csv("df_count.csv", index=False)
    print("finished")


if __name__ == '__main__':
    main()
