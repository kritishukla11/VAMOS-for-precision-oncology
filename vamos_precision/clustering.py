import logging
logger = logging.getLogger(__name__)

def run_clustering(**kwargs):
    """
    Auto-wrapped from 2_Density_Based_Clustering.py.
    kwargs can pass configuration and file paths.
    """
    #filenames
    input_data = "alphafold_scores_above70.csv"
    output_data = "ready_for_ML_processed.csv"

    #import all necessary packages
    import pandas as pd
    import numpy as np
    from sklearn.neighbors import NearestNeighbors
    import matplotlib.pyplot as plt
    from sklearn.cluster import DBSCAN
    from collections import Counter
    import seaborn as sns
    import statistics as st
    from sklearn.metrics.pairwise import pairwise_distances
    import sklearn
    import itertools
    import os
    from kneed import KneeLocator, DataGenerator as dg
    import numpy
    import warnings

    #turn off warnings
    warnings.filterwarnings('ignore')

    #read in data
    df1 = pd.read_csv(input_data)

    #first do one gene, here VPS13D
    df2 = df1[df1['gene_name']=='VPS13D']
    df = pd.DataFrame().assign(x =df2['x_coord'], y =df2['y_coord'], z =df2['z_coord'])

    nbrs = NearestNeighbors(n_neighbors=5).fit(df)
    neigh_dist, neigh_ind = nbrs.kneighbors(df)
    sort_neigh_dist = np.sort(neigh_dist, axis=0)
    k_dist = sort_neigh_dist[:, 4]
    k_dist_list = k_dist.tolist()
    length = list(range(len(k_dist_list)))


    x, y = length, k_dist_list
    kl = KneeLocator(x, y, curve="convex")
    epsilon = int(k_dist_list[kl.knee])/2

    clusters = DBSCAN(eps= epsilon, min_samples=6).fit(df)
    cluster_list = clusters.labels_.tolist()
    df.insert(0, 'cluster', cluster_list)
    df['gene'] = 'VPS13d'
    df = df.reset_index()

    cluster_list2 = [*set(cluster_list)]

    # do one cluster

    df_cluster_loop = df[df['cluster'] == 0]
    dist_matrix = pd.DataFrame().assign(x=df_cluster_loop['x'], y=df_cluster_loop['y'], z=df_cluster_loop['z'])
    pw_dist = pairwise_distances(dist_matrix)
    avg_dist = numpy.average(pw_dist)
    df_cluster_loop["dist"] = avg_dist
    df_cluster_loop['num_members'] = len(df_cluster_loop['gene'])
    df_cluster = df_cluster_loop.copy()

    #do the rest
    cluster_list2.remove(0)

    for i in cluster_list2:
        df_cluster_loop = df[df['cluster'] == i]
        dist_matrix = pd.DataFrame().assign(x=df_cluster_loop['x'], y=df_cluster_loop['y'], z=df_cluster_loop['z'])
        pw_dist = pairwise_distances(dist_matrix)
        avg_dist = numpy.average(pw_dist)
        df_cluster_loop["dist"] = avg_dist
        df_cluster_loop['num_members'] = len(df_cluster_loop['gene'])

        dfs = [df_cluster, df_cluster_loop]
        df_cluster = pd.concat(dfs)

    #get all genes that have 3 or more mutations
    genes = df1['gene_name'].to_list()
    genes.remove('VPS13D')

    res  = []
    for x in set(genes):
        if genes.count(x) >= 3:
            res.append(x)

    #do all genes 
    #possible error: small number of genes may not have a clear "knee", epsilon for those has to be picked manually
    for i in res:
        subset = df1[df1["gene_name"] == i]
        xyz = pd.DataFrame().assign(x =subset['x_coord'], y =subset['y_coord'], z =subset['z_coord'])
        nbrs = NearestNeighbors(n_neighbors=5).fit(xyz)
        neigh_dist, neigh_ind = nbrs.kneighbors(xyz)
        sort_neigh_dist = np.sort(neigh_dist, axis=0)
        k_dist = sort_neigh_dist[:, 4]
        k_dist_list = k_dist.tolist()
        length = list(range(len(k_dist_list)))

        x, y = length, k_dist_list
        kl = KneeLocator(x, y, curve="convex")
        epsilon = int(k_dist_list[kl.knee])/2

        clusters = DBSCAN(eps= epsilon, min_samples=6).fit(xyz)
        cluster_list = clusters.labels_.tolist()
        xyz.insert(0, 'cluster', cluster_list)
        xyz['gene'] = i
        xyz = xyz.reset_index()

        cluster_list2 = [*set(cluster_list)]

        for x in cluster_list2:
            df_cluster_loop = xyz[xyz['cluster'] == x]
            dist_matrix = pd.DataFrame().assign(x=df_cluster_loop['x'], y=df_cluster_loop['y'], z=df_cluster_loop['z'])
            pw_dist = pairwise_distances(dist_matrix)
            avg_dist = numpy.average(pw_dist)
            df_cluster_loop["dist"] = avg_dist
            df_cluster_loop['num_members'] = len(df_cluster_loop['gene'])

            dfs = [df_cluster, df_cluster_loop]
            df_cluster = pd.concat(dfs)

    #merge both dfs to get all info
    df_final = pd.merge(df1,df_cluster,on=['x','y','z','gene','index1'])

    #save
    df_final.to_csv(output_data)
    logger.info('Completed run_clustering.')
    return True
