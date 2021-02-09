"""
Main package funcion of pancovcorrect
"""

# standard import
import inspect
import os
import sys

# dependencies import
import snakemake

# project import
from . import __name__ as package_name
from . import cli


def main(**config):
    """
    Main function of pancovcorrect

    Parameter:

    config: A dict pass to snakemake as config.
    Value associate to 'snakemake' is parse with snakemake cli parser.

    Return:

    True if workflow execution was successful.
    """
    snakemake_parser = snakemake.get_argument_parser()
    params = vars(snakemake_parser.parse_args(config["snakemake"]))
    del config["snakemake"]

    # Found Snakefile
    params["snakefile"] = os.path.join(
        os.path.dirname(__file__), f"{package_name}.snk"
    )

    # Add parameter with config
    params["config"] = config

    # Add target
    params["targets"] = ["all"]

    # Clean not use argument
    args_set = set(inspect.signature(snakemake.snakemake).parameters.keys())
    params = {k: v for k, v in params.items() if k in args_set}

    return snakemake.snakemake(**params)


if __name__ == "__main__":
    # Get parameter in assign failed stop
    if args := cli.get_argument(sys.argv[1:]):
        sys.exit(int(not main(**args)))
    else:
        sys.exit(1)
