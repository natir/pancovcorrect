"""
Main package funcion of pancovcorrect
"""

# standard import
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
    params = config["snakemake"]
    del config["snakemake"]

    # Found Snakefile
    params.extend(
        [
            "--snakefile",
            os.path.join(os.path.dirname(__file__), f"{package_name}.snk"),
        ]
    )

    # Add parameter with config
    params.append("--config")
    for (name, value) in config.items():
        params.append(f"{name}={value}")

    # Add target
    params.extend(["--", "all"])

    return snakemake.main(params)


if __name__ == "__main__":
    # Get parameter in assign failed stop
    if args := cli.get_argument(sys.argv[1:]):
        sys.exit(int(not main(**args)))
    else:
        sys.exit(1)
