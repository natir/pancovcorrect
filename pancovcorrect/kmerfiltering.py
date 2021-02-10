"""
Function use to perform kmer filtering
"""

# standard import
import csv
import itertools
from collections import Counter


def kmerfiltering(input_path, output_path, **config):
    """
    Filter kmer
    """
    kmer_counts = read_kmc_dump(input_path)

    with open(output_path, "w") as out:
        if "reference" in input_path:
            for (kmer, _) in generate_kmc_dump(input_path["reference"]):
                print(f">references\n{kmer}", file=out)

        for kmer in valid_kmer(kmer_counts, config):
            print(f">reads\n{kmer}", file=out)


def valid_kmer(kmer_counts, config):
    """
    Generate only valid kmer
    """
    for forward in kmer_counts.keys():
        reverse = __rev_comp(forward)
        tt_counts = kmer_counts[forward] + kmer_counts[reverse]

        if tt_counts > 0:
            ratio = min(
                kmer_counts[forward] / tt_counts,
                kmer_counts[reverse] / tt_counts,
            )
        else:
            ratio = 0

        if ratio < config["forward_min_ratio"]:
            continue

        if tt_counts < config["min_abundance"]:
            continue

        if tt_counts > config["max_abundance"]:
            yield forward
        else:
            preds = [
                kmer_counts[pred_mer] + kmer_counts[__rev_comp(pred_mer)]
                for pred_mer in around_kmer(forward, -1)
                if pred_mer in kmer_counts
                or __rev_comp(pred_mer) in kmer_counts
            ]
            succs = [
                kmer_counts[succ_mer] + kmer_counts[__rev_comp(succ_mer)]
                for succ_mer in around_kmer(forward, 1)
                if succ_mer in kmer_counts
                or __rev_comp(succ_mer) in kmer_counts
            ]

            if preds and succs:
                coverage = min(max(preds), max(succs))
            else:
                coverage = 0

            if tt_counts > (coverage * config["guess_coverage_ratio"]):
                yield forward


def generate_kmc_dump(path):
    """
    Generate (kmer, count) from kmc dump file
    """
    with open(str(path)) as inp:
        reader = csv.reader(inp, delimiter=" ")
        for row in reader:
            yield (row[0], int(row[1]))


def read_kmc_dump(path):
    """
    Convert kmc dump file in Counter
    """
    data = Counter()

    for (kmer, count) in generate_kmc_dump(path):
        data[kmer] = count

    return data


def around_kmer(kmer, pos):
    """
    Generate kmer around actual kmer
    """
    if pos < 0:
        suffix = kmer[:pos]
        for prefix in __generate_all_seq(abs(pos)):
            new_kmer = prefix + suffix
            if new_kmer != kmer:
                yield new_kmer

    elif pos > 0:
        prefix = kmer[pos:]
        for suffix in __generate_all_seq(pos):
            new_kmer = prefix + suffix
            if new_kmer != kmer:
                yield new_kmer

    else:
        yield kmer


def __generate_all_seq(length):
    """
    Generate all sequence with length
    """
    for sequence in itertools.product(["A", "C", "T", "G"], repeat=length):
        yield "".join(sequence)


def __rev_comp(seq):
    """
    Return reverse complement of sequence
    """

    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(reversed([complement[n] for n in seq]))
