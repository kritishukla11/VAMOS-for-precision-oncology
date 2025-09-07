import logging
logger = logging.getLogger(__name__)

def compute_priority_scores(**kwargs):
    """
    Auto-wrapped from 7_Priority_scores.py.
    kwargs can pass configuration and file paths.
    """
    #name files
    cell_RNA = 'WILLIAMS_ESR1_TARGETS_DN-BreastOnly-logodds.csv'
    cell_CRISPR = 'ESR1_CRISPR-BreastOnly-logodds.csv')
    tumor_RNA = 'TCGA-GOZGIT_ESR1_TARGETS_UP-BreastOnly-logodds.csv'
    scores = 'ESR1_DN_priority_score.csv'

    #import packages
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np

    #read in files
    df1 = pd.read_csv(cell_RNA)
    df2 = pd.read_csv(cell_CRISPR)
    df3 = pd.read_csv(tumor_RNA)

    #get priority score
    list1 = df1['id'].tolist()
    list2 = df2['id'].tolist()
    list3 = df3['id'].tolist()
    clust = [*set(list1 + list2 + list3)]

    k_list = []

    for i in clust:
        k = 0
        n = 0
        if i in list1:
            k = k + df1[df1['id']==i].reset_index()['log_odds'][0]
            n = n+1
        if i in list2:
            k = k + df2[df2['id']==i].reset_index()['log_odds'][0]
            n = n+1
        if i in list3:
            k = k + df3[df3['id']==i].reset_index()['log_odds'][0]
            n = n+1

        k_list.append(k/n)

    #make dataframe
    df = pd.DataFrame({'cluster': clust, 'priority': k_list}).sort_values(by='priority',ascending=False)

    #filter to keep only high scores, remove outlier clusters, keep highest score per gene
    df = df[df['priority']>6]
    df = df[df['clust']!='-1'].reset_index(drop=True)
    df[['gene','clust']] = df['cluster'].str.split('_',expand=True).reset_index(drop=True)
    df = df.loc[df.groupby("gene")["priority"].idxmax()].sort_values(by='priority',ascending=False).reset_index()

    #save
    df.to_csv(scores)
    logger.info('Completed compute_priority_scores.')
    return True
