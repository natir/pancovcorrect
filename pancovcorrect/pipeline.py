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

    if "id" in wcs.keys():
        return f"{config['working_dir']}/reads/{k}/{wcs.id}.fasta"

    return f"{config['working_dir']}/reads/{k}/{wcs.filename}.fasta"


def select_kmer_input(config):
    """
    Generate input for kmer selection step with and without reference
    """

    if config["reference"] != "None":
        return {
            "reads": f"{config['working_dir']}/kmersets/{{k}}/{{id}}.tsv",
            "reference": f"{config['working_dir']}/kmersets/{{k}}/reference.tsv",
        }

    return {"reads": f"{config['working_dir']}/kmersets/{{k}}/{{id}}.tsv"}
