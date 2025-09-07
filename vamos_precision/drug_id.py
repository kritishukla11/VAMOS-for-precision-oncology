import logging
logger = logging.getLogger(__name__)

def identify_drugs(**kwargs):
    """
    Auto-wrapped from 8_Drug_ID.py.
    kwargs can pass configuration and file paths.
    """
    #files
    drug_data = 'Drug_sensitivity_AUC_(CTD^2)_subsetted.csv'
    struc_data = 'WILLIAMS_ESR1_TARGETS_DN-afterML.csv'
    final_drugs = "sig_drugs_PIK3CA_1.csv"

    #cluster of interest
    clust = 'PIK3CA_1'

    #import packages
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import ranksums
    import warnings

    #turn off warnings
    warnings.filterwarnings('ignore')

    #read in data
    auc_og_esr1 = pd.read_csv(drug_data)
    auc_og_esr1.columns = auc_og_esr1.columns.str.replace(r"\s*\(.*?\)", "", regex=True)

    #get cell lines in and not in cluster
    df = pd.read_csv(struc_data)
    df['id'] = df['gene'] + '_' + df['cluster'].astype(str)
    sub_in = df[df['id']==clust]
    sub_out = df[df['id']!=clust]

    list_mut = [*set(sub_in['DepMap_ID'].tolist())]
    list_norm1 = [*set(sub_out['DepMap_ID'].tolist())]
    list_norm = [i for i in list_norm1 if i not in list_mut]

    auc_mut = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_mut)]
    auc_norm = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_norm)]

    #get all drugs with info for the cluster
    drugs = [*set(auc_mut_stacked.dropna().index.tolist())]

    #get p values for drugs with more of an effect
    drug_list = []
    p_value_list = []

    for i in drugs:
        if auc_mut_stacked.dropna().loc[i]['Values'].mean() < auc_norm_stacked.dropna().loc[i]['Values'].mean():
            statistic, p_value = ranksums(auc_mut_stacked.loc[i].dropna()['Values'],auc_norm_stacked.loc[i].dropna()['Values'])
            if p_value < 0.05:
                drug_list.append(i)
                p_value_list.append(p_value)
                mean_list.append(auc_mut_stacked.loc[i].dropna()['Values'].mean())
        else:
            continue

    #make dataframe       
    sig_drugs = pd.DataFrame({'drug': drug_list, 'p_value': p_value_list}).sort_values(by='p_value').reset_index(drop=True)

    #save
    sig_drugs.to_csv(final_drugs, index=False)
    logger.info('Completed identify_drugs.')
    return True
