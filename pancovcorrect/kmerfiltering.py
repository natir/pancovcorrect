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
        for kmer in valid_kmer(kmer_counts, config):
            print("f>1\n{kmer}", file=out)


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
        elif tt_counts < config["min_abundance"]:
            continue
        elif tt_counts > config["max_abundance"]:
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


def read_kmc_dump(path):
    """
    Read kmc dump
    """
    data = Counter()

    with open(str(path)) as fh:
        reader = csv.reader(fh, delimiter=" ")
        for row in reader:
            data[row[0]] = int(row[1])

    return data


def contain_homopolymer(kmer, length=4):
    """
    Return true if kmer contains homopolymer
    """
    currentRepeatLength = 0
    lastBase = "N"
    for base in kmer:
        if base == lastBase:
            currentRepeatLength += 1
            if currentRepeatLength == length:
                return True
        else:
            currentRepeatLength = 0
        lastBase = base
    return False


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
    for p in itertools.product(["A", "C", "T", "G"], repeat=length):
        yield "".join(p)


def __rev_comp(seq):
    """
    Return reverse complement of sequence
    """

    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(reversed([complement[n] for n in seq]))
