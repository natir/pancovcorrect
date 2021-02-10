"""
Function use full for snakemake pipeline
"""


def uncorrected_path(wcs, config):
    """
    Generate uncorrected path from wildcard and config
    """

    if wcs.k == str(config["kmer_sizes"][0]):
        k = "raw"
    else:
        k = config["kmer_sizes"][config["kmer_sizes"].index(wcs.k) - 1]

    return f"{config['working_dir']}/reads/{k}/{wcs.filename}.fasta"


def select_kmer_input(wcs, config):
    """
    Generate input for kmer selection step with and without reference
    """

    inputs = dict()

    inputs["reads"] = f"{config['working_dir']}/kmersets/{{k}}/{{filename}}.tsv"
    if config['reference'] != 'None':
        inputs["reference"] = f"{config['working_dir']}/kmersets/{{k}}/reference.tsv"

    return inputs
