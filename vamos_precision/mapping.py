import logging
logger = logging.getLogger(__name__)

def run_mapping(**kwargs):
    """
    Auto-wrapped from 1_Mapping_to_Alphafold.py.
    kwargs can pass configuration and file paths.
    """
    #file names
    CCLE_data = "CCLE_mutations.csv"
    uniprot_mapping_data = "uniprot_to_gene.csv"
    alphafold_3D = "alphafold_3Dcoord_mapped.csv"
    alphafold_struc = "alphafold_struc_mapped.csv"
    alphafold_scores = "alphafold_scores_mapped.csv"

    #function for Alphafold filename automation
    def bar(s):
        return 'AF-' + s + '-F1-model_v3.cif'

    #importing all necessary packages
    #this is for mapping
    from unipressed import IdMappingClient
    import time
    #pandas
    import pandas as pd
    #cif reading
    import sys
    from gemmi import cif
    #plotting
    import matplotlib.pyplot as plt
    #joining lists, opening files
    import itertools
    import os
    #warnings
    import warnings

    #turn off warnings
    warnings.filterwarnings('ignore')

    #import CCLE data
    depmap = pd.read_csv(CCLE_data)

    #import Uniprot ID to gene name mapping file
    info = pd.read_csv(uniprot_mapping_data)

    #map in Uniprot IDs
    depmap['uniprot'] = "Unknown"

    for i in range(len(depmap)):
        if depmap['gene'][i] in info.index:
            depmap['uniprot'][i] = info.loc[depmap['gene'][i]][0]

    #do one one gene first, here VPS13D
    #subset VPS data
    VPS = metadata[metadata['gene_name'] == 'VPS13D']
    #map to Alphafold 
    doc = cif.read('AF-Q5THJ4-F1-model_v3.cif')
    block = doc.sole_block()

    atom_id = block.find_loop('_atom_site.id')
    x_coord = block.find_loop('_atom_site.Cartn_x')
    y_coord = block.find_loop('_atom_site.Cartn_y')
    z_coord = block.find_loop('_atom_site.Cartn_z')

    atom_id = [eval(i) for i in atom_id]
    atom_id = pd.to_numeric(atom_id)

    x_coord = [eval(i) for i in x_coord] 
    y_coord = [eval(i) for i in y_coord] 
    z_coord = [eval(i) for i in z_coord] 

    xyz_coord_dic = {'atom_mut':atom_id,'x_coord':x_coord, 'y_coord':y_coord, 'z_coord':z_coord}
    xyz_coord_df = pd.DataFrame(xyz_coord_dic)
    #merge xyz with gene info
    tot = pd.merge(xyz_coord_df, VPS, on='atom_mut')

    #get all genes
    gene_list = [*set(depmap['gene_name'].to_list())]

    #map to Alphafold
    no_struct_info = []

    for x in gene_list:
        y = depmap[depmap['gene'] == x].reset_index(drop=True)['uniprot_id'][0]

        if os.path.isfile(bar(y)):
            with open(bar(y), 'r') as f1:

                doc = cif.read(bar(y))
                block = doc.sole_block()
                atom_id = block.find_loop('_atom_site.id')
                x_coord = block.find_loop('_atom_site.Cartn_x')
                y_coord = block.find_loop('_atom_site.Cartn_y')
                z_coord = block.find_loop('_atom_site.Cartn_z')

                atom_id = [eval(i) for i in atom_id]
                atom_id = pd.to_numeric(atom_id)
                x_coord = [eval(i) for i in x_coord] 
                y_coord = [eval(i) for i in y_coord] 
                z_coord = [eval(i) for i in z_coord] 

                xyz_coord_dic = {'atom_mut':atom_id,'x_coord':x_coord, 'y_coord':y_coord, 'z_coord':z_coord}
                xyz_coord_df = pd.DataFrame(xyz_coord_dic)

                df_new = pd.merge(xyz_coord_df, sub, on='atom_mut')

                dfs = [tot, df_new]

                tot = pd.concat(dfs)

        else:
            no_struct_info.append(x)
            continue

    tot = tot.reset_index(drop=True)

    #save all mappings
    tot.to_csv(alphafold_3D)

    #add atom type info

    #again first for VPS13D
    doc = cif.read('AF-Q5THJ4-F1-model_v3.cif')
    block = doc.sole_block()
    pos1 = block.find_loop('_struct_conf.beg_auth_seq_id')
    struc1 = block.find_loop('_struct_conf.conf_type_id')
    pos = [eval(i) for i in pos1]
    pos = pd.to_numeric(pos)
    struc = [i for i in struc1]
    struc_dic = {'position': pos,'structure':struc}
    struc_df = pd.DataFrame(struc_dic)
    pos2 = block.find_loop('_struct_conf.end_auth_seq_id')
    struc2 = block.find_loop('_struct_conf.id')
    pos_ = [eval(i) for i in pos2]
    pos_ = pd.to_numeric(pos_)
    struc_ = [i for i in struc2]
    struc2_dic = {'position': pos_,'structure':struc_}
    struc2_df = pd.DataFrame(struc2_dic)
    dfs = [struc_df,struc2_df]
    df = pd.concat(dfs)
    df = df.reset_index(drop = True)
    df['structure_pos'] = df['structure'].str.replace('([A-Z]+)', '')
    df['structure_type'] = df['structure'].str.extract('([A-Z]+)')
    df = df.drop(columns=['structure_pos'])
    df = df.drop(columns=['structure'])
    df['disorder'] = df['structure_type'].map(lambda x: x == "STRN")
    df['position'] = df['position'].astype(int)
    df['position'] = df['position'].astype(int)

    #add to 3D coordinate df
    sub = tot[tot["gene_name"] == 'VPS13D']
    merged = pd.merge(df, sub, on=['position'])

    #now for all genes
    for i in genes:
        sub = tot[tot['gene_name'] == i]
        sub = sub.reset_index(drop=True)

        y = sub['uniprot_id'][0]

        if os.path.isfile(bar(y)):
            with open(bar(y), 'r') as f1:

                doc = cif.read(bar(y))
                block = doc.sole_block()
                pos1 = block.find_loop('_struct_conf.beg_auth_seq_id')
                struc1 = block.find_loop('_struct_conf.conf_type_id')
                pos = [eval(i) for i in pos1]
                pos = pd.to_numeric(pos)
                struc = [i for i in struc1]
                struc_dic = {'position': pos,'structure':struc}
                struc_df = pd.DataFrame(struc_dic)

                pos2 = block.find_loop('_struct_conf.end_auth_seq_id')
                struc2 = block.find_loop('_struct_conf.id')
                pos_ = [eval(i) for i in pos2]
                pos_ = pd.to_numeric(pos_)
                struc_ = [i for i in struc2]
                struc2_dic = {'position': pos_,'structure':struc_}
                struc2_df = pd.DataFrame(struc2_dic)

                dfs = [struc_df,struc2_df]
                structure_df = pd.concat(dfs)

                structure_df['structure_pos'] = structure_df['structure'].str.replace('([A-Z]+)', '')
                structure_df['structure_type'] = structure_df['structure'].str.extract('([A-Z]+)')
                structure_df = structure_df.drop(columns=['structure_pos'])
                structure_df = structure_df.drop(columns=['structure'])
                structure_df['disorder'] = structure_df['structure_type'].map(lambda x: x == "STRN")
                structure_df['position'] = structure_df['position'].astype(int)

                merged2 = pd.merge(structure_df, sub, on=['position'])

                dfs = [merged2,merged]
                merged = pd.concat(dfs)

    #save this
    merged.to_csv(alphafold_struct)

    #add confidence score info
    #formatting dataframe to add confidence scores by isolating protein residue
    xyz[['firstcol', 'mut_residue', 'thirdcol']] = xyz['Protein_Change'].str.extract('([A-Za-z]+)(\d+\.?\d*)([A-Za-z]*)', expand = True)
    xyz = xyz.drop(columns=['firstcol','thirdcol'])
    xyz = xyz.dropna(subset=['mut_residue'])
    xyz = xyz.reset_index(drop=True)

    #add scores

    scores_list = []

    for i in range(len(xyz)):
        x = merged['uniprot_id'][i] 
        if os.path.isfile(bar(x)):
            with open(bar(x), 'r') as f1:
                doc = cif.read(bar(x))
                block = doc.sole_block()
                residue = block.find_loop('_ma_qa_metric_local.label_seq_id')
                score = block.find_loop('_ma_qa_metric_local.metric_value')
                scorelist = list(score)
                residuelist = list(residue)
                residuelist = [eval(i) for i in residuelist]
                scorelist = list(score)
                y = int(xyz['mut_residue'][i])
                if y in residuelist:
                    scores_list.append(scorelist[y-1])
                else:
                    scores_list.append('no') 

    xyz['score'] = scores_list

    #save this
    xyz.to_csv(alphafold_scores)
    logger.info('Completed run_mapping.')
    return True
