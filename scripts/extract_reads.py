from subprocess import run
import argparse
import struct
import mmap
import sys
import os


class FastqIndexElement:
    def __init__(self, sequence_offset, quality_offset, length):
        self.sequence_offset = sequence_offset
        self.quality_offset = quality_offset
        self.length = length

    def __str__(self):
        return str(self.sequence_offset) + "\t" + str(self.quality_offset) + "\t" + str(self.length)


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


'''
faidx format:
-------------
NAME	Name of this reference sequence
LENGTH	Total length of this reference sequence, in bases
OFFSET	Offset in the FASTA/FASTQ file of this sequence's first base
LINEBASES	The number of bases on each line
LINEWIDTH	The number of bytes in each line, including the newline
QUALOFFSET	Offset of sequence's first quality within the FASTQ file
'''
def load_fastq_index(faidx_path):
    name_to_offset = dict()
    index_elements = list()

    with open(faidx_path, 'r') as file:
        for line in file:
            data = line.strip().split('\t')

            name = data[0]
            length = int(data[1])
            sequence_offset = int(data[2])
            quality_offset = int(data[5])

            index_element = FastqIndexElement(sequence_offset, quality_offset, length)

            index_elements.append(index_element)

            name_to_offset[name] = len(index_elements) - 1

    return name_to_offset, index_elements


def load_query_ids(query_ids_path):
    ids = list()
    with open(query_ids_path, 'r') as file:
        for line in file:
            ids.append(int(line.strip()))

    return ids


def extract_bytes_from_file(mmap_file_object, offset, n_bytes):
    s = None

    data = mmap_file_object[offset:offset+n_bytes]

    s = bytes(data)

    return s


# Use a system call to samtools faidx to build the index
def build_index(path):
    index_path = path + ".fai"
    if os.path.exists(index_path):
        sys.stderr.write("Index exists\n")
        if os.path.getmtime(index_path) < os.path.getmtime(path):
            sys.stderr.write("WARNING: fasta/q index path is older than file itself, may be out of date: " + index_path + "\n")
    else:
        sys.stderr.write("No index found, indexing... ")
        arguments = ["samtools", "faidx", path]
        run(arguments, check=True)
        sys.stderr.write("Done\n")

    sys.stderr.flush()

    return index_path


def main(csv_path, fastq_path, query_ids_path):
    faidx_path = fastq_path + ".fai"

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
