from modules.Fastx import *

from subprocess import run
import argparse
import random
import struct
import mmap
import sys
import os


def main(fastq_path, n_reads, min_length):
    faidx_path = build_index(fastq_path)

    name_to_offset, index_elements = load_fastq_index(faidx_path=faidx_path)
    items = list(name_to_offset.items())

    random.shuffle(items)

    output_path = os.path.splitext(fastq_path)[0] + "_random_" + str(n_reads) + "_min-length-" + str(min_length) + ".fastq"
    sys.stderr.write("Writing to: " + output_path + '\n')

    i = 0
    n = 0
    with open(fastq_path, 'rb') as input_file, open(output_path, 'wb') as output_file:
        mm = mmap.mmap(input_file.fileno(), 0, prot=mmap.PROT_READ)

        while n < n_reads and i < len(items):
            name,offset_index = items[i]

            index_element = index_elements[offset_index]

            if index_element.length >= min_length:
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

                # Only increment n if the read passes the length requirement
                n += 1

            # Always move on to next sequence
            i += 1

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
    parser.add_argument(
        "--min_length",
        type=int,
        required=False,
        default=0,
        help="skip any reads less than this length"
    )

    args = parser.parse_args()

    main(
        fastq_path=args.fastq,
        n_reads=args.n,
        min_length=args.min_length,
    )
