from modules.IterativeHistogram import IterativeHistogram

from matplotlib import pyplot
import math
import sys
import os


def get_chromosome_ordering(name):
    ordinal = 0
    ordinal_string = name.split("chr")[-1]

    has_alpha = False
    has_numeric = False
    for c in ordinal_string:
        if c.isalpha():
            has_alpha = True
        elif c.isnumeric():
            has_numeric = True

    if has_alpha and not has_numeric:
        weight = 1
        for c in ordinal_string:
            ordinal += ord(c) * weight
            weight *= 0.1

    elif has_numeric and not has_alpha:
        ordinal = int(ordinal_string)

    else:
        exit("ERROR: region name is both numeric and alphabetical")

    return ordinal


"""
 Parse a PAF file: https://github.com/lh3/miniasm/blob/master/PAF.md using the following data:
 
 Col (1-based)       Data                        Type
 -------------       ----                        ----
 1                   "Query sequence name"       string (assumed to be numeric for this project, using ONT reads)
 6                   "Target sequence name"      string
 7                   "Target sequence length"    int
 8                   "Target start..."           int
 9                   "Target end..."             int
 12                  "Mapping quality"           int
"""
def main():
    paf_path = "/home/ryan/data/nanopore/human/test/cross_strand_analysis/HG002_fc2_palindromes_run1_VS_chm13.paf"

    mapping_distributions_per_chromosome = dict()
    max_region_size = 0
    bin_size = 500_000

    with open(paf_path, 'r') as file:
        for l,line in enumerate(file):
            data = line.strip().split()

            region = data[5]
            region_length = int(data[6])
            start = int(data[7])
            stop = int(data[8])

            midpoint = (stop + start) / 2

            if region not in mapping_distributions_per_chromosome:
                n_bins = int(math.ceil(float(region_length)/bin_size))

                mapping_distributions_per_chromosome[region] = \
                    IterativeHistogram(start=0, stop=region_length, n_bins=n_bins)

                if region_length > max_region_size:
                    max_region_size = region_length

                get_chromosome_ordering(region)

            mapping_distributions_per_chromosome[region].update(midpoint)

    n_rows = int(math.ceil(float(len(mapping_distributions_per_chromosome))/2))
    figure, axes = pyplot.subplots(nrows=n_rows, ncols=2, sharey=True)

    for i,item in enumerate(sorted(mapping_distributions_per_chromosome.items(), key=lambda x: get_chromosome_ordering(x[0]))):
        region,iterative_histogram = item

        row_index = int(math.floor(float(i)/n_rows))
        column_index = i % n_rows

        print(row_index, column_index)

        frequencies = iterative_histogram.get_histogram()
        bin_centers = iterative_histogram.get_bin_centers()
        x = [float(i)/2 for i in range(len(frequencies))]

        axes[column_index][row_index].bar(x=x, height=frequencies)
        axes[column_index][row_index].set_ylabel(region, rotation=0, labelpad=30)
        axes[column_index][row_index].set_xlim([0,250])

    axes[column_index][-1].set_xlabel("Coordinate (Mbp)")

    # pyplot.subplots_adjust(hspace=1.5)

    figure.set_size_inches(18,24)

    # pyplot.show()
    # pyplot.close()

    pyplot.savefig("mapping_distribution.png", dpi=200)


if __name__ == "__main__":
    main()
