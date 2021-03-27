from matplotlib import pyplot
import argparse
import numpy
import os


def parse_read_lengths_from_file(path):
    lengths = list()

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            line = line.strip()
            lengths.append(int(line))

    return lengths


def main(a_path, b_path):
    a = parse_read_lengths_from_file(a_path)
    b = parse_read_lengths_from_file(b_path)

    a_frequencies,bins = numpy.histogram(a, bins=500, range=(0,100000), normed=True)
    b_frequencies,bins = numpy.histogram(b, bins=500, range=(0,100000), normed=True)

    centers = (bins[:-1] + bins[1:]) / 2
    width = centers[1]-centers[0]

    axes = pyplot.axes()
    l1 = axes.plot(centers, a_frequencies, label=os.path.basename(a_path))
    l2 = axes.plot(centers, b_frequencies, label=os.path.basename(b_path))
    axes.legend()

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-a",
        type=str,
        required=True,
        help="path of file containing a list of read lengths, one on each line"
    )
    parser.add_argument(
        "-b",
        type=str,
        required=True,
        help="path of file containing a list of read lengths, one on each line"
    )

    args = parser.parse_args()

    main(
        a_path=args.a,
        b_path=args.b
    )
