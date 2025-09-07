import logging
logger = logging.getLogger(__name__)

def compute_logodds(**kwargs):
    """
    Auto-wrapped from 6_Logodds.py.
    kwargs can pass configuration and file paths.
    """
    #read_data
    input_file = "GOZGIT_ESR1_TARGETS_UP-afterML.csv"
    output_file = "WILLIAMS_ESR1_TARGETS_DN-BreastOnly-logodds.csv"

    #import necessary packages
    import pandas as pd 
    import numpy as np

    #read in file
    df = pd.read_csv(input_file)

    #genelist
    genes = [*set(df['gene'].to_list())]

    #calculate log_odds
    gene_list = []
    cluster_list = []
    log_odds_list = []

    for x in genes:
        sub = df[df['gene']==x]
        N = len(sub)

        array_clust = sub['cluster'].to_numpy()
        array_class1 = sub['class'].to_numpy()

        joint_counts = pd.crosstab(array_clust,array_class1, rownames=['cluster'], colnames=['class'])
        joint_prob = joint_counts/N

        sub_cluster_list = [*set(sub['cluster'].to_list())]

        for i in sub_cluster_list:
            sub_cluster_df = sub[sub['cluster']==i]

            if len([*set(sub_cluster_df['label'].to_list())]) >= 3:

                if 0 in set(sub['class']) and 1 in set(sub_cluster_df['class']) and 0 in set(sub_cluster_df['class']):
                    d = sum(joint_counts[0])-joint_counts[0][i]
                    b = joint_counts[0][i]
                    c = sum(joint_counts[1])-joint_counts[1][i]
                    a = joint_counts[1][i]

                    odds = np.true_divide(a,c+0.00000001)
                    OR = np.true_divide(np.multiply(a,d),np.multiply(b,c))
                    LOR = np.log(OR)

                    gene_list.append(x)
                    cluster_list.append(i)
                    log_odds_list.append(LOR)
                elif 0 in set(sub['class']) and 1 in set(sub_cluster_df['class']):
                    d = sum(joint_counts[0])
                    b = 0
                    c = sum(joint_counts[1])-joint_counts[1][i]
                    a = joint_counts[1][i]

                    odds = np.true_divide(a,c)
                    OR = np.true_divide(np.multiply(a,d),(np.multiply(b,c)+0.001))
                    LOR = np.log(OR)

                    gene_list.append(x)
                    cluster_list.append(i)
                    log_odds_list.append(LOR)

                elif 1 in set(sub_cluster_df['class']):
                    d = 0
                    b = 0
                    c = sum(joint_counts[1])-joint_counts[1][i]
                    a = joint_counts[1][i]

                    odds = np.true_divide(a,c)
                    OR = np.true_divide(np.multiply(a,d),(np.multiply(b,c)+0.001))
                    LOR = np.log(OR)

                    gene_list.append(x)
                    cluster_list.append(i)
                    log_odds_list.append(LOR)

                else:
                    continue

            else:
                continue

    #make dataframe
    df2= pd.DataFrame(list(zip(gene_list, cluster_list,log_odds_list)),
                      columns =['gene', 'cluster','log_odds'])

    #get rid of outliers
    df3 = df2[df2['log_odds']<100]
    df3 = df3[df3['log_odds']>-100]
    df3 = df3.dropna()

    #save
    df3.to_csv(output_file)
    logger.info('Completed compute_logodds.')
    return True
