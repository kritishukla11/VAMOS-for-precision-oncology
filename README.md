# VAMOS_precision_oncology

Developed by Kriti Shukla, kritis@unc.edu

# What this does:
Using this workflow, you can obtain whole proteome structural information from Alphafold, create dense variant clusters, identify which clusters are associated with a given pathway, run statistical analysis on these clusters to identify likelihood of association, and find drugs differentially associated with those clusters in cell lines.

# Instructions to run:
Data : All input data is publicly available on DepMap, TCGA, CTRP, Alphafold, and MSigDB.

Environment: Set up conda environment using the provided YAML file

# Known issues:
This workflow uses the kneed package for automatic finding of the epsilon parameter during density based clustering. If your input data does not have a clear "elbow" or "knee," this automatic process will fail and a manual epsilon parameter will need to be assigned for that protein. 

# Citations:
Please cite the following paper: 

## Installation

    python -m venv .venv
    source .venv/bin/activate
    pip install -e .

## CLI Usage

    vamos-precision --help
    vamos-precision -c config.yaml map
    vamos-precision -c config.yaml cluster
    vamos-precision -c config.yaml gsea-prep
    vamos-precision -c config.yaml gsea-run   # requires R
    vamos-precision -c config.yaml ml
    vamos-precision -c config.yaml logodds
    vamos-precision -c config.yaml prioritize
    vamos-precision -c config.yaml drug-id
    vamos-precision -c config.yaml permute

## Examples
See `examples/quickstart.ipynb` and toy CSVs in `examples/` to get started.
