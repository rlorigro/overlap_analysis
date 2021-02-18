from modules.Fastq import *

from subprocess import run
import argparse
import struct
import mmap
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
                id_to_name.extend([None]*(len(id_to_name) - id + 1))

            id_to_name[id] = name

    return id_to_name


def load_query_ids(query_ids_path):
    ids = list()
    with open(query_ids_path, 'r') as file:
        for line in file:
            ids.append(int(line.strip()))

    return ids


def main(csv_path, fastq_path, query_ids_path):
    faidx_path = build_index(fastq_path)

    id_to_name = generate_id_to_name_mapping(csv_path=csv_path)

    name_to_offset, index_elements = load_fastq_index(faidx_path=faidx_path)

    queries = load_query_ids(query_ids_path=query_ids_path)

    output_path = os.path.splitext(fastq_path)[0] + "_" + os.path.splitext(os.path.basename(query_ids_path))[0] + ".fastq"
    sys.stderr.write("Writing to: " + output_path + '\n')

    with open(fastq_path, 'rb') as input_file, open(output_path, 'wb') as output_file:
        mm = mmap.mmap(input_file.fileno(), 0, prot=mmap.PROT_READ)

        for id in queries:
            if id >= len(id_to_name):
                sys.stderr.write("WARNING: skipping id " + str(id) + " because it is greater than max id in shasta ReadSummary.csv\n")
                continue

            name = id_to_name[id]
            offset_index = name_to_offset[name]
            index_element = index_elements[offset_index]
            print(name, index_element)

            s = extract_bytes_from_file(mmap_file_object=mm,
                                        offset=index_element.sequence_offset,
                                        n_bytes=index_element.length)

            q = extract_bytes_from_file(mmap_file_object=mm,
                                        offset=index_element.quality_offset,
                                        n_bytes=index_element.length)

            output_file.write(b'@')
            output_file.write(name.encode('utf-8'))
            output_file.write(b'\n')
            output_file.write(s)
            output_file.write(b'\n')
            output_file.write(b'+')
            output_file.write(b'\n')
            output_file.write(q)
            output_file.write(b'\n')

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fastq",
        type=str,
        required=True,
        help="path of file containing FASTA/FASTQ sequence"
    )
    parser.add_argument(
        "--summary_csv",
        type=str,
        required=True,
        help="path of ReadSummary.csv containing mapping from readId (shasta) to readName (fastq)"
    )
    parser.add_argument(
        "--ids",
        type=str,
        required=True,
        help="path of file containing 1 id per line to be queried"
    )

    args = parser.parse_args()

    main(
        fastq_path=args.fastq,
        csv_path=args.summary_csv,
        query_ids_path=args.ids,
    )
