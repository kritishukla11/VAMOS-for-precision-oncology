import argparse
from .utils import setup_logging
from .config import Config
from .mapping import run_mapping
from .clustering import run_clustering
from .gsea import prep_gsea_inputs, run_gsea_r_wrapper
from .ml import run_ml
from .logodds import compute_logodds
from .prioritization import compute_priority_scores
from .drug_id import identify_drugs
from .permutation import run_permutation_analysis

def main():
    parser = argparse.ArgumentParser(prog="vamos-precision", description="VAMOS Precision Oncology Pipeline")
    parser.add_argument("-v", "--verbose", action="count", default=1, help="Increase verbosity (-v, -vv)")
    parser.add_argument("-c", "--config", help="Path to YAML configuration", default=None)

    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("map", help="Mapping to AlphaFold")
    sub.add_parser("cluster", help="Density-based clustering")
    sub.add_parser("gsea-prep", help="Prepare GSEA inputs (Python)")
    sub.add_parser("gsea-run", help="Run GSEA (R)")
    sub.add_parser("ml", help="Machine learning pipeline")
    sub.add_parser("logodds", help="Compute log-odds metrics")
    sub.add_parser("prioritize", help="Compute priority scores")
    sub.add_parser("drug-id", help="Identify drug associations")
    sub.add_parser("permute", help="Permutation analysis")

    args, unknown = parser.parse_known_args()
    setup_logging(args.verbose)
    cfg = Config(args.config) if args.config else Config()

    if args.cmd == "map":
        run_mapping(**cfg.get("mapping", {}))
    elif args.cmd == "cluster":
        run_clustering(**cfg.get("clustering", {}))
    elif args.cmd == "gsea-prep":
        prep_gsea_inputs(**cfg.get("gsea", {}))
    elif args.cmd == "gsea-run":
        run_gsea_r_wrapper(**cfg.get("gsea", {}))
    elif args.cmd == "ml":
        run_ml(**cfg.get("ml", {}))
    elif args.cmd == "logodds":
        compute_logodds(**cfg.get("logodds", {}))
    elif args.cmd == "prioritize":
        compute_priority_scores(**cfg.get("prioritization", {}))
    elif args.cmd == "drug-id":
        identify_drugs(**cfg.get("drug_id", {}))
    elif args.cmd == "permute":
        run_permutation_analysis(**cfg.get("permutation", {}))

if __name__ == "__main__":
    main()
