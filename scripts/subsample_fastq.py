from modules.Fastq import *

from subprocess import run
import argparse
import random
import struct
import mmap
import sys
import os


def main(fastq_path, n_reads):
    faidx_path = build_index(fastq_path)

    name_to_offset, index_elements = load_fastq_index(faidx_path=faidx_path)
    items = list(name_to_offset.items())

    random.shuffle(items)

    output_path = os.path.splitext(fastq_path)[0] + "_random.fastq"
    sys.stderr.write("Writing to: " + output_path + '\n')

    with open(fastq_path, 'rb') as input_file, open(output_path, 'wb') as output_file:
        mm = mmap.mmap(input_file.fileno(), 0, prot=mmap.PROT_READ)

        for i in range(n_reads):
            name,offset_index = items[i]

            index_element = index_elements[offset_index]

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
        "--n",
        type=int,
        required=True,
        help="number of reads to subsample"
    )

    args = parser.parse_args()

    main(
        fastq_path=args.fastq,
        n_reads=args.n,
    )
