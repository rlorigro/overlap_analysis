from modules.Fastq import *
from modules.Kmer import *
import mmap


def main():
    fastq_path = "/home/ryan/data/nanopore/human/test/cross_strand_analysis/GM24385_2_Guppy_3.6.0_prom_random_reads.fastq"

    faidx_path = build_index(fastq_path)

    output_directory = os.path.join(os.path.dirname(fastq_path), "plots")

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    name_to_offset, index_elements = load_fastq_index(faidx_path=faidx_path)

    k = 10
    density = 0.1

    with open(fastq_path, 'r') as file:
        mm = mmap.mmap(file.fileno(), 0, prot=mmap.PROT_READ)

        for name,i in name_to_offset.items():
            index_element = index_elements[i]

            if index_element.length < 1000:
                continue

            s = extract_bytes_from_file(mmap_file_object=mm,
                                        offset=index_element.sequence_offset,
                                        n_bytes=index_element.length)

            marker_iterator = MarkerIterator(s, k=k, density=density)

            n_markers = 0
            while marker_iterator.next():
                # print(marker_iterator.kmer_string(), marker_iterator.kmer_number(), marker_iterator.marker_number())
                n_markers += 1

            print(len(s), n_markers, float(len(s))*density)


if __name__ == "__main__":
    main()
