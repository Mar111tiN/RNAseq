import sys

# logging
sys.stderr = open(s.log[0], "w")

import pandas as pd


def get_matrix(s):

    counts = [
        pd.read_table(
            f, index_col=0, usecols=[0, 1], header=None, skiprows=4
        ) for f in s.input
    ]

    for t, sample in zip(counts, s.params.samples):
        t.columns = [sample]

    matrix = pd.concat(counts, axis=1)
    matrix.index.name = "gene"
    # collapse technical replicates
    matrix = matrix.groupby(matrix.columns, axis=1).sum()
    matrix.to_csv(s.output[0], sep="\t")
