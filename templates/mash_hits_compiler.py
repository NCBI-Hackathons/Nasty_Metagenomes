#!/usr/bin/env python3

"""
Purpose
-------

This module intends to retrieve the AMR fasta sequences from mash.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``sample_id`` : Sample Identification string.
    - e.g.: ``'SampleA'``
- ``amr_reference`` : A fasta file path.
    - e.g.: ``'AMR_CDS.fasta'``
- ``mash_results`` : mash hits
    - e.g.: ``'.screen'``

Generated output
----------------

-  A fasta file per contig (given the minimum contig size
"""

__version__ = "0.0.1"
__build__ = "27082019"
__process__ = "mash_hits_compiler-nf"

import os
from itertools import groupby
import csv


if __file__.endswith(".command.sh"):
    SAMPLE_ID = '$sample_id'
    REFERENCE = '$amr_reference'
    MASH_RESULTS = '$mash_results'
    print("Running {} with parameters:".format(
        os.path.basename(__file__)))
    print("SAMPLE_ID: {}".format(SAMPLE_ID))
    print("ASSEMBLY: {}".format(REFERENCE))
    print("MIN_SIZE: {}".format(MASH_RESULTS))


def main(sample_id, reference, mash_results):
    """Main executor of the split_fasta template.

    Parameters
    ----------
    sample_id : str
        Sample Identification string.
    assembly : list
        Assembly file.
    min_size : int
        Minimum contig size."""

    print("Starting script")

    amr_hits = []

    with open(mash_results, "r") as mash:
        data = mash.readlines()
        for line in data:
            amr_hits.append(line.split('\\t')[4])

    print(len(amr_hits))

    f_open = open(reference, "r")
    results_file = open(sample_id + '_amr.fasta', 'w')

    entry = (x[1] for x in groupby(f_open, lambda line: line[0] == ">"))

    for header in entry:

        header_str = header.__next__()[1:].strip()
        seq = "".join(s.strip() for s in entry.__next__())

        print(header_str.split(' ')[0])

        if header_str.split(' ')[0] in amr_hits:
            results_file.write(">" + header_str + "\\n" + seq + "\\n")

    f_open.close()
    results_file.close()


if __name__ == '__main__':
    main(SAMPLE_ID, REFERENCE, MASH_RESULTS)
