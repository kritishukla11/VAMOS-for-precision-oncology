# VAMOS_precision_oncology

Developed by Kriti Shukla, kritis@unc.edu

# What this does:
Using this workflow, you can apply the previously published VAMOS ML model for drug discovery in precision oncology.

# Instructions to run:
Data : All input data is publicly available on DepMap, TCGA, CTRP, Alphafold, and MSigDB.

Environment: Set up conda environment using the provided YAML file

# Known issues:
This workflow uses the kneed package for automatic finding of the epsilon parameter during density based clustering. If your input data does not have a clear "elbow" or "knee," this automatic process will fail and a manual epsilon parameter will need to be assigned for that protein. 

# Citations:
Please cite the following BioRXiv link (now accepted at ASPCT CTS): https://www.biorxiv.org/content/10.1101/2025.03.14.643357v1 

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
