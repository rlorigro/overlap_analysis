import argparse
import sys
import os


def generate_id_to_name_mapping(csv_path):
    id_to_name = list()

    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split(',')

            id = int(data[0])
            name = data[1]

            if len(id_to_name) <= id:
                id_to_name.extend([None]*(abs(len(id_to_name) - id) + 1))

            id_to_name[id] = name

    return id_to_name


def convert(overlaps_csv_path, id_to_name):
    output_path = os.path.splitext(overlaps_csv_path)[0] + "_named.csv"

    with open(overlaps_csv_path, 'r') as file, open(output_path, 'w') as out_file:
        for l,line in enumerate(file):
            if l == 0:
                out_file.write(line.replace("Id", "Name"))
                continue

            data = line.strip().split(',')

            id1 = int(data[0])
            id2 = int(data[1])

            converted_line = id_to_name[id1] + ',' + id_to_name[id2] + ',' + data[2] + '\n'

            out_file.write(converted_line)


def main(overlaps_csv_path, summary_csv_path):
    id_to_name = generate_id_to_name_mapping(summary_csv_path)
    convert(overlaps_csv_path, id_to_name)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--overlaps",
        type=str,
        required=True,
        help="path of file containing FASTA/FASTQ sequence"
    )
    parser.add_argument(
        "--summary",
        type=str,
        required=True,
        help="path of ReadSummary.csv containing mapping from readId (shasta) to readName (fastq)"
    )

    args = parser.parse_args()

    main(
        overlaps_csv_path=args.overlaps,
        summary_csv_path=args.summary,
    )
