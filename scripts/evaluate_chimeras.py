from subprocess import run
import argparse
import sys
import os


def run_evaluation(reference_path, fastq_path, s3_path, dry, n_threads):
    ref_name = "_".join(os.path.basename(reference_path).split('.')[:-1])
    fastq_name = "_".join(os.path.basename(fastq_path).split('.')[:-1])

    output_name = fastq_name + "_VS_" + ref_name
    output_dir = os.path.abspath(output_name)
    print("\nSTARTING: " + output_name)

    if os.path.exists(output_name):
        print("Output already exists. Skipping: " + output_name)
        return
    else:
        os.mkdir(output_dir)

    project_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    build_dir = os.path.join(project_path, "build")
    executable_path = os.path.join(build_dir, "filter_chimeras_from_paf")

    if not os.path.exists(executable_path):
        exit("ERROR: executable path not found in build directory of project")

    print(executable_path)

    if n_threads is None:
        n_threads = os.cpu_count()

    minimap_args = [
        "minimap2",
        "-x", "map-ont",
        "--secondary=no",
        "-n", "10",
        "-K", "10g",
        "-k", "17",
        "-t", str(n_threads),
        reference_path,
        fastq_path
    ]

    print(' '.join(minimap_args))
    minimap_output_path = os.path.join(output_dir, output_name + ".paf")

    if not dry:
        with open(minimap_output_path, 'w') as file:
            run(minimap_args, check=True, stdout=file, stderr=sys.stderr, universal_newlines=True)

    chimera_exe_args = [
        executable_path,
        "--paf_path", minimap_output_path
    ]

    if not dry:
        run(chimera_exe_args, check=True)


def main(reference_path, fastq_paths, s3_path, dry, n_threads):
    for path in fastq_paths.split(','):
        run_evaluation(reference_path=reference_path, fastq_path=path, s3_path=s3_path, dry=dry, n_threads=n_threads)
        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ref","-r",
        type=str,
        required=True,
        help="Path of reference fasta to use for alignment"
    )

    parser.add_argument(
        "--fastq","-i",
        type=str,
        required=True,
        help="comma-separated list of fastq file paths to be aligned"
    )

    parser.add_argument(
        "--s3",
        type=str,
        default=None,
        required=False,
        help="AWS s3 URI of directory to upload to, if not specified, no attempt to upload is made"
    )

    parser.add_argument(
        "--dry",
        dest="dry",
        action="store_true",
        required=False,
        help="echo the commands instead of executing them"
    )

    parser.add_argument(
        "--n_threads",
        type=int,
        default=None,
        required=False,
        help="echo the commands instead of executing them"
    )

    args = parser.parse_args()

    main(
        reference_path=args.ref,
        fastq_paths=args.fastq,
        s3_path=args.s3,
        dry=args.dry,
        n_threads=args.n_threads
    )
