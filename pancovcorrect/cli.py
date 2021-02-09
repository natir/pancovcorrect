"""
Parse command line argument for pancovcorrect
"""

# standard import
import argparse
import logging

# project import
from . import __name__ as package_name


def get_argument(args):
    """
    Use argparse to analyse args and return a dict with argument
    """

    parser = argparse.ArgumentParser(
        prog=package_name,
        description="Pancovcorrect long reads with iterative kmer length",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Required params
    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        required=True,
        help="List of read inputs file",
    )
    parser.add_argument(
        "-o",
        "--outputs",
        nargs="+",
        required=True,
        help="List of read outputs file",
    )
    parser.add_argument(
        "-w",
        "--working-dir",
        type=str,
        required=True,
        help="A directory where working file is write",
    )

    parser.add_argument(
        "-k",
        "--kmer-sizes",
        nargs="+",
        type=int,
        default=[11, 13, 15],
        help="List of kmer size used",
    )

    parser.add_argument(
        "-r",
        "--reference",
        type=str,
        help="Path to reference genome (optional)",
    )

    # Kmer filter params
    parser.add_argument(
        "-m",
        "--min-abundance",
        type=int,
        default=10,
        help="Kmer with abundance lower than this parameter is discard",
    )
    parser.add_argument(
        "-M",
        "--max-abundance",
        type=int,
        default=50,
        help="Kmer with abundance upper than this parameter is keep",
    )
    parser.add_argument(
        "-f",
        "--forward-min-ratio",
        type=float,
        default=0.3,
        help="Kmer with forward / total abundance lower than this parameter is discard",
    )
    parser.add_argument(
        "-c",
        "--guess-coverage-ratio",
        type=float,
        default=0.7,
        help="Kmer abundance between min and max abundance with coverage upper than guess coverage times this parameter is keep",
    )

    # Snakemake params
    parser.add_argument(
        "-t",
        "--max-threads",
        type=int,
        default=48,
        help="Max threads per rules",
    )

    parser.add_argument(
        "snakemake",
        nargs=argparse.REMAINDER,
        help="all argument after '--' is pass to snakemake",
    )

    args = parser.parse_args(args)

    # Edit input
    args.kmer_sizes.sort()

    args.snakemake = [val for val in args.snakemake if val != "--"]
    # Check input
    if len(args.inputs) != len(args.outputs):
        logging.error("inputs and outputs need have same number of values")
        return None

    if len(args.kmer_sizes) != len(set(args.kmer_sizes)):
        logging.error("each kmer size must be uniq")
        return None

    if args.min_abundance >= args.max_abundance:
        logging.error("min-abundance need be lower than max-abundance")
        return None

    if args.forward_min_ratio < 0 or args.forward_min_ratio > 1:
        logging.error("forward-min-ratio need be in intervale [0 and 1]")
        return None

    if args.guess_coverage_ratio < 0 or args.guess_coverage_ratio > 1:
        logging.error("guess-coverage-ratio need be in intervale [0 and 1]")
        return None

    return vars(args)
