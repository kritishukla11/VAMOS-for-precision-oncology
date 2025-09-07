import logging
logger = logging.getLogger(__name__)

def run_permutation_analysis(**kwargs):
    """
    Auto-wrapped from 9_Permutation_Analysis.py.
    kwargs can pass configuration and file paths.
    """
    #files
    drug_data = 'Drug_sensitivity_AUC_(CTD^2)_subsetted.csv'
    struc_data = 'WILLIAMS_ESR1_TARGETS_DN-afterML.csv'
    output_file = 'PIK3CA_permutations_df.csv'

    #import packages
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy.stats import ranksums
    import random
    import numpy as np
    import statistics
    import warnings

    #turn off warnings
    warnings.filterwarnings('ignore')

    #read in data
    auc_og = pd.read_csv(drug_data)
    auc_og.columns = auc_og_esr1.columns.str.replace(r"\s*\(.*?\)", "", regex=True)

    #load in structural data 
    df = pd.read_csv(struc_data)
    df['id'] = df['gene'] + '_' + df['cluster'].astype(str)

    #drugs of interest PIK3CA UP
    drug_list = ['ISOX:BORTEZOMIB','MYRICETIN','BRD8899','SMER-3','MK 2206','PX-12','17-AAG','WAY-362450','PARTHENOLIDE',
            'ML083','TANESPIMYCIN:GEMCITABINE','PRIMA-1']

    #get info for a cluster, repeat this for all clusters
    p_value_list = []

    for drug in drug_list:
        index_list = ['Unnamed: 0',drug]
        auc_og_esr1 = auc_og[index_list]

        clust = 'PIK3CA_1'

        sub4 = df3[df3['id']==clust]
        sub5 = df3[df3['id']!=clust]

        list_mut = [*set(sub4['DepMap_ID'].tolist())]
        list_norm1 = [*set(sub5['DepMap_ID'].tolist())]
        list_norm = [i for i in list_norm1 if i not in list_mut]

        auc_mut = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_mut)]
        auc_norm = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_norm)]

        auc_mut_stacked = auc_mut.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
        auc_mut_stacked['mut'] = 0.8

        auc_norm_stacked = auc_norm.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
        auc_norm_stacked['mut']=0.2

        auc = pd.concat([auc_norm_stacked.dropna(),auc_mut_stacked.dropna()])

        statistic, p_value = ranksums(auc_mut_stacked.dropna()['Values'],auc_norm_stacked.dropna()['Values'])

        p_value_list.append(p_value)


    #make dataframe
    PIK3CA_df = pd.DataFrame({'drug': drug_list, 'VAMOS_clust': p_value_list})

    # permutations 100
    mean_list = []

    for drug in drug_list:
        index_list = ['Unnamed: 0',drug]
        auc_og_esr1 = auc_og[index_list]

        p_value_list = []

        for i in range(0,100):
            random_numbers = random.sample(list(df3[df3['gene']=='PIK3CA'].index), 4)

            sub4 = df3.loc[random_numbers]
            sub5 = df3.drop(index=random_numbers)

            list_mut = [*set(sub4['DepMap_ID'].tolist())]
            list_norm1 = [*set(sub5['DepMap_ID'].tolist())]
            list_norm = [i for i in list_norm1 if i not in list_mut]

            auc_mut = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_mut)]
            auc_norm = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_norm)]

            auc_mut_stacked = auc_mut.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
            auc_mut_stacked['mut'] = 0.8

            auc_norm_stacked = auc_norm.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
            auc_norm_stacked['mut']=0.2

            auc = pd.concat([auc_norm_stacked.dropna(),auc_mut_stacked.dropna()])

            statistic, p_value = ranksums(auc_mut_stacked.dropna()['Values'],auc_norm_stacked.dropna()['Values'])

            p_value_list.append(p_value)

        cleaned_list = [x for x in p_value_list if not pd.isna(x)]
        mean_list.append(statistics.mean(cleaned_list))

    #add to dataframe
    PIK3CA_df['random_100_permutations'] = mean_list

    # permutations 1000
    mean_list = []

    for drug in drug_list:
        index_list = ['Unnamed: 0',drug]
        auc_og_esr1 = auc_og[index_list]

        p_value_list = []

        for i in range(0,1000):
            random_numbers = random.sample(list(df3[df3['gene']=='PIK3CA'].index), 4)

            sub4 = df3.loc[random_numbers]
            sub5 = df3.drop(index=random_numbers)

            list_mut = [*set(sub4['DepMap_ID'].tolist())]
            list_norm1 = [*set(sub5['DepMap_ID'].tolist())]
            list_norm = [i for i in list_norm1 if i not in list_mut]

            auc_mut = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_mut)]
            auc_norm = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_norm)]

            auc_mut_stacked = auc_mut.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
            auc_mut_stacked['mut'] = 0.8

            auc_norm_stacked = auc_norm.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
            auc_norm_stacked['mut']=0.2

            auc = pd.concat([auc_norm_stacked.dropna(),auc_mut_stacked.dropna()])

            statistic, p_value = ranksums(auc_mut_stacked.dropna()['Values'],auc_norm_stacked.dropna()['Values'])

            p_value_list.append(p_value)

        cleaned_list = [x for x in p_value_list if not pd.isna(x)]
        mean_list.append(statistics.mean(cleaned_list))

    #add to dataframe
    PIK3CA_df['random_1000_permutations'] = mean_list

    # permutations 10000
    mean_list = []

    for drug in drug_list:
        index_list = ['Unnamed: 0',drug]
        auc_og_esr1 = auc_og[index_list]

        p_value_list = []

        for i in range(0,10000):
            random_numbers = random.sample(list(df3[df3['gene']=='PIK3CA'].index), 4)

            sub4 = df3.loc[random_numbers]
            sub5 = df3.drop(index=random_numbers)

            list_mut = [*set(sub4['DepMap_ID'].tolist())]
            list_norm1 = [*set(sub5['DepMap_ID'].tolist())]
            list_norm = [i for i in list_norm1 if i not in list_mut]

            auc_mut = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_mut)]
            auc_norm = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_norm)]

            auc_mut_stacked = auc_mut.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
            auc_mut_stacked['mut'] = 0.8

            auc_norm_stacked = auc_norm.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
            auc_norm_stacked['mut']=0.2

            auc = pd.concat([auc_norm_stacked.dropna(),auc_mut_stacked.dropna()])

            statistic, p_value = ranksums(auc_mut_stacked.dropna()['Values'],auc_norm_stacked.dropna()['Values'])

            p_value_list.append(p_value)

        cleaned_list = [x for x in p_value_list if not pd.isna(x)]
        mean_list.append(statistics.mean(cleaned_list))


    #add to dataframe
    PIK3CA_df['random_10000_permutations'] = mean_list


    # permutations 100000
    mean_list = []

    for drug in drug_list:
        index_list = ['Unnamed: 0',drug]
        auc_og_esr1 = auc_og[index_list]

        p_value_list = []

        for i in range(0,100000):
            random_numbers = random.sample(list(df3[df3['gene']=='PIK3CA'].index), 4)

            sub4 = df3.loc[random_numbers]
            sub5 = df3.drop(index=random_numbers)

            list_mut = [*set(sub4['DepMap_ID'].tolist())]
            list_norm1 = [*set(sub5['DepMap_ID'].tolist())]
            list_norm = [i for i in list_norm1 if i not in list_mut]

            auc_mut = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_mut)]
            auc_norm = auc_og_esr1[auc_og_esr1['Unnamed: 0'].isin(list_norm)]

            auc_mut_stacked = auc_mut.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
            auc_mut_stacked['mut'] = 0.8

            auc_norm_stacked = auc_norm.set_index('Unnamed: 0').T.melt(ignore_index=False,var_name="Original_Column", value_name="Values")
            auc_norm_stacked['mut']=0.2

            auc = pd.concat([auc_norm_stacked.dropna(),auc_mut_stacked.dropna()])

            statistic, p_value = ranksums(auc_mut_stacked.dropna()['Values'],auc_norm_stacked.dropna()['Values'])

            p_value_list.append(p_value)

        cleaned_list = [x for x in p_value_list if not pd.isna(x)]
        mean_list.append(statistics.mean(cleaned_list))

    #add to dataframe
    PIK3CA_df['random_100000_permutations'] = mean_list


    #save
    PIK3CA_df.to_csv(output_file,index=False)
    logger.info('Completed run_permutation_analysis.')
    return True
