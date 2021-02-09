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
