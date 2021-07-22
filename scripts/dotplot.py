from modules.Fastx import *

from multiprocessing import cpu_count
from matplotlib import pyplot
from matplotlib import cm
from subprocess import run
import argparse
import numpy
import mmap
import os


def align_minimap(ref_sequence_path, reads_sequence_path, max_threads=None, output_dir=None, preset="map-ont", k=15, all_chains=False, additional_args=None):
    """
    Given a reference file and reads file align using minimap, generating a
    :param ref_sequence_path:
    :param reads_sequence_path:
    :param output_dir:
    :return:
    """

    if max_threads is None:
        max_threads = max(1, cpu_count() - 2)

    max_threads = str(max_threads)

    ref_sequence_path = os.path.abspath(ref_sequence_path)
    reads_sequence_path = os.path.abspath(reads_sequence_path)

    print(reads_sequence_path, ref_sequence_path)

    ref_sequence_filename_prefix = os.path.basename(ref_sequence_path)
    ref_sequence_filename_prefix = os.path.splitext(ref_sequence_filename_prefix)[0]

    input_filename_prefix = os.path.basename(reads_sequence_path)
    input_filename_prefix = os.path.splitext(input_filename_prefix)[0]
    output_filename_prefix = input_filename_prefix + "_VS_" + ref_sequence_filename_prefix

    # ---- Minimap -----------

    output_filename = output_filename_prefix + ".sam"
    output_file_path = os.path.join(output_dir, output_filename)

    arguments = ["minimap2", "-a", "-t", max_threads, "-x", preset, "-k", str(k)]

    if all_chains:
        arguments.append("-P")

    if additional_args is not None:
        arguments += additional_args

    arguments += [ref_sequence_path, reads_sequence_path]

    print("\nRUNNING: ", " ".join(arguments))

    with open(output_file_path, "w") as output_file:
        print("REDIRECTING TO: ", output_file_path, "\n")
        run(arguments, cwd=output_dir, stdout=output_file, check=True)

    # ---- Sort SAM ----------

    input_filename = output_filename
    output_filename = output_filename_prefix + ".sorted.sam"

    arguments = ["samtools", "sort", input_filename, "-@", max_threads, "-O", "SAM", "-o", output_filename]

    print("\nRUNNING: ", " ".join(arguments))

    run(arguments, cwd=output_dir, check=True)

    os.remove(os.path.join(output_dir, input_filename))

    return os.path.join(output_dir, output_filename)


def parse_cigar_as_tuples(cigar_string):
    operations = list()

    length_string = ""

    for i,c in enumerate(cigar_string):
        if c.isalpha():
            operations.append((c, int(length_string)))
            length_string = ""
        if c.isnumeric():
            length_string += c

    return operations


def is_reference_move(cigar_type):
    if cigar_type == 'M':
        return True
    elif cigar_type == 'I':
        return False
    elif cigar_type == 'D':
        return True
    elif cigar_type == '=':
        return True
    elif cigar_type == 'X':
        return True
    elif cigar_type == 'S':
        return False
    elif cigar_type == 'H':
        return False


def is_query_move(cigar_type):
    if cigar_type == 'M':
        return True
    elif cigar_type == 'I':
        return True
    elif cigar_type == 'D':
        return False
    elif cigar_type == '=':
        return True
    elif cigar_type == 'X':
        return True
    elif cigar_type == 'S':
        return True
    elif cigar_type == 'H':
        return False


def get_ref_alignment_length(cigar_operations):
    l = 0
    for item in cigar_operations:
        l += item[1]*is_reference_move(item[0])

    return l


def generate_dot_plot(sam_path, sequence_length, title, output_directory):
    figure = pyplot.figure()
    axes = pyplot.axes()

    viridis = cm.get_cmap('viridis', 256)

    colors = {'M':"blue",
              'I':"orange",
              'D':"orange",
              '=':"blue",
              'X':"purple",
              'S':"black",
              'H':"black"}

    with open(sam_path, 'r') as file:
        for l,line in enumerate(file):
            if line[0] == "@":
                continue

            data = line.strip().split()

            name = data[0]
            flags = int(data[1])
            start_index = int(data[3])
            cigar_string = data[5]
            quality = int(data[4]) if data[4] != "*" else 0

            # Some self-alignments have an insane number of secondaries, cant plot them without crashing
            if l > 20 and quality < 1:
                continue

            is_reverse = (flags & 0x10 == 0x10)

            print(name, is_reverse, start_index, cigar_string[:min(10,len(cigar_string))])

            cigar_operations = parse_cigar_as_tuples(cigar_string)

            ref_index = start_index - 1
            query_index = 0

            if is_reverse:
                ref_index += get_ref_alignment_length(cigar_operations)
                query_index += sequence_length
                cigar_operations = reversed(cigar_operations)

            for operation in cigar_operations:
                x1 = ref_index
                x2 = x1 + int(is_reference_move(operation[0]))*(1-2*int(is_reverse))*operation[1]
                y1 = query_index
                y2 = y1 + int(is_query_move(operation[0]))*(1-2*int(is_reverse))*operation[1]

                # print(operation, is_reference_move(operation[0]), is_query_move(operation[0]), "(%d,%d) -> (%d,%d)" % (x1,y1,x2,y2))

                ref_index = x2
                query_index = y2

                if operation[0] != "S" and operation[0] != "H":
                    axes.plot([x1,x2], [y1,y2], color=colors[operation[0]], linewidth=0.5)

    axes.set_aspect('equal')
    axes.set_ylim([0,sequence_length])
    axes.set_xlim([0,sequence_length])
    axes.set_title(title)

    pyplot.savefig(os.path.join(output_directory, title+"_self_alignment.png"), dpi=200)

    # pyplot.show()
    # pyplot.close()


def generate_quality_plot(qualities_string, title, output_directory):
    figure = pyplot.figure()
    axes = pyplot.axes()

    window_size = 100
    step_size = 6.0/float(window_size)

    x_kernel = numpy.arange(start=-3,stop=3+step_size,step=step_size)
    kernel = (2.71828**(-(x_kernel**2)/2))/numpy.sqrt(2*3.14159)
    kernel /= numpy.sum(kernel)

    qualities = numpy.array([(10**(float(q-33)/-10)) for q in map(ord, list(qualities_string))], dtype=numpy.float)

    smoothed_qualities = numpy.convolve(qualities, kernel, 'valid')

    axes.plot(smoothed_qualities)

    axes.set_ylim([0,0.7])

    figure.set_size_inches(12,4)
    axes.set_title(title)

    pyplot.savefig(os.path.join(output_directory, title+"_quality.png"), dpi=200)


def main(fastq_path):
    faidx_path = build_index(fastq_path)

    output_directory = os.path.splitext(fastq_path)[0] + "_plots/"

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    name_to_offset, index_elements = load_fastq_index(faidx_path=faidx_path)

    with open(fastq_path, 'r') as file:
        mm = mmap.mmap(file.fileno(), 0, prot=mmap.PROT_READ)

        for name,i in name_to_offset.items():
            index_element = index_elements[i]

            if index_element.length < 1000:
                continue

            s = extract_bytes_from_file(mmap_file_object=mm,
                                        offset=index_element.sequence_offset,
                                        n_bytes=index_element.length)

            q = extract_bytes_from_file(mmap_file_object=mm,
                                        offset=index_element.quality_offset,
                                        n_bytes=index_element.length)

            output_file_path = os.path.join(output_directory, name + ".fastq")

            with open(output_file_path, 'wb') as output_file:
                output_file.write(b'@')
                output_file.write(name.encode('utf-8'))
                output_file.write(b'\n')
                output_file.write(s)
                output_file.write(b'\n')
                output_file.write(b'+')
                output_file.write(b'\n')
                output_file.write(q)
                output_file.write(b'\n')

            sam_path = align_minimap(
                output_file_path,
                output_file_path,
                k=11,
                all_chains=True,
                additional_args=["--rev-only"],
                output_dir=output_directory,
                )

            generate_dot_plot(sam_path, sequence_length=len(s), title=name, output_directory=output_directory)

            os.remove(output_file_path)
            os.remove(sam_path)

            generate_quality_plot(q.decode("utf-8"), title=name, output_directory=output_directory)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        type=str,
        required=True,
        help="path of fastq file"
    )

    args = parser.parse_args()

    main(
        fastq_path=args.i,
    )
