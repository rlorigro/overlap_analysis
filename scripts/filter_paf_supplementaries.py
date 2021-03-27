import sys


def filter_by_matches(paf_path):
    """
    PAF has columns (0-based):
        0 = read_name
        9 = n_matches
    """

    alignments_per_read = dict()

    with open(paf_path, 'r') as file:
        for l,line in enumerate(file):
            tokens = line.strip().split()

            read_name = tokens[0]

            if read_name in alignments_per_read:
                prev_tokens = alignments_per_read[read_name].strip().split()
                n_matches = int(tokens[9])
                prev_n_matches = int(prev_tokens[9])

                if n_matches > prev_n_matches:
                    alignments_per_read[read_name] = line
            else:
                alignments_per_read[read_name] = line

    output_path = ''.join(paf_path.strip().split('.')[:-1]) + "_filtered-most-matches.paf"

    with open(output_path, 'w') as file:
        for line in alignments_per_read.values():
            file.write(line)

    return


def main(paf_path):
    filter_by_matches(paf_path=paf_path)

    return


if __name__ == "__main__":

    if len(sys.argv) != 2:
        exit("ERROR: need to provide 1 argument: path of input PAF")

    main(sys.argv[1])
