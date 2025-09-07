import logging
logger = logging.getLogger(__name__)

def prep_gsea_inputs(**kwargs):
    """
    Auto-wrapped from 3_GSEA_python_prepfiles.py.
    kwargs can pass configuration and file paths.
    """
    #input file names
    tsvfile = 'GOZGIT_ESR1_TARGETS_UP.tsv'
    columnname = 'GOZGIT_ESR1_TARGETS_UP'
    gct_file_name = 'GOZGIT_ESR1_TARGETS_UP_GE.csv'
    gmtfile = 'GOZGIT_ESR1_TARGETS_UP'
    genesetname = 'GOZGIT_ESR1_TARGETS_UP.csv'
    gene_expression_data = "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
    sample_info = "sample_info.csv"

    #importing necessary packages
    import pandas as pd
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    import warnings

    #turn off warnings
    warnings.filterwarnings('ignore')

    #load in our data
    samples = pd.read_csv(sample_info)
    GE = pd.read_csv(gene_expression_data)

    #reformatting the gene expression dataframe to match the format we want
    GE.rename(columns={'Unnamed: 0':'DepMap_ID'}, inplace=True)
    GE = GE.merge(samples[['DepMap_ID','CCLE_Name']],on='DepMap_ID')
    cols = GE.columns

    new_set = []
    new_set.append(cols[0])
    for i in cols[1:-1]:
        new_set.append(i.split(' ')[0])

    new_set.append(cols[-1])

    GE.columns = new_set

    #load in the gene set
    gene_set = pd.read_table(tsvfile)

    #get the list of genes we care about
    gene_set_genes = [i for i in gene_set.loc[16,columnname].split(',') if (len(i) > 0) and (i in GE.columns)]
    gene_set_genes_noGE = [i for i in gene_set.loc[16,columnname].split(',') if (len(i) > 0) and (i not in GE.columns)]

    #subset our gene expression dataframe to only the genes we care about
    gene_set_GE = GE[['CCLE_Name']+gene_set_genes]

    #reformat the gene expression file again and check to see it looks good
    gct_file = gene_set_GE
    gct_file = gct_file.set_index('CCLE_Name').transpose().reset_index()
    gct_file['Name'] = gene_set_genes
    gct_file['Description'] = gene_set_genes
    column_to_move = gct_file.pop("Description")
    gct_file.insert(0, "Description", column_to_move)
    column_to_move = gct_file.pop("Name")
    gct_file.insert(0, "Name", column_to_move)
    gct_file = gct_file.drop(columns=["index"])
    gct_file = gct_file.set_index('Name')
    gct_file = gct_file.drop(columns=['Description'])

    #save the file to computer
    gct_file.to_csv(gct_file_name)

    #load in the gene set file
    gene_set_df = pd.read_csv(gmtfile,sep='\t',header=None).T

    #make the top row into the column name (reformatting file)
    new_header = gene_set_df.iloc[0] #grab the first row for the header
    gene_set_df = gene_set_df[1:] #take the data less the header row
    gene_set_df.columns = new_header #set the header row as the df header

    #grab only the rows we care about
    gene_set_df = gene_set_df.iloc[1:len(gene_set_df)-1].reset_index(drop=True)

    #save the file to computer
    gene_set_df.to_csv(genesetname,index=False)
    logger.info('Completed prep_gsea_inputs.')
    return True


from .utils import run_r_script, ensure_dir
from pathlib import Path

def run_gsea_r_wrapper(output_dir: str | Path = "results/gsea", args: list[str] | None = None):
    """Run the bundled R GSEA script via Rscript. Ensure R is installed."""
    ensure_dir(output_dir)
    script_path = Path(__file__).parent / "r" / "GSEA_final.R"
    return run_r_script(script_path, args or [])
