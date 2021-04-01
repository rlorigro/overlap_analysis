from subprocess import run
import argparse
import sys
import os


def run_evaluation(shasta_path, fastq_path, s3_path, dry):
    fastq_path = os.path.abspath(fastq_path)
    fastq_filename = os.path.basename(fastq_path)
    fastq_prefix = os.path.splitext(fastq_filename)[0]
    fastq_directory = os.path.dirname(fastq_path)
    output_name = fastq_prefix + "_palindromes/"
    output_directory = os.path.join(fastq_directory, output_name)

    shasta_args = [
        shasta_path,
        "--input", fastq_path,
        "--assemblyDirectory", output_directory,
        "--Reads.minReadLength", "10000",
        "--command", "filterReads",
        "--Reads.palindromicReads.detectOnFastqLoad",
        "--Reads.palindromicReads.nearDiagonalFractionThreshold", "0",
        "--Reads.palindromicReads.maxSkip", "150",
        "--Reads.palindromicReads.maxDrift", "150",
    ]

    if dry:
        shasta_args = ["echo"] + shasta_args

    run(shasta_args, check=True)

    aws_args = [
        "aws", "s3", "cp",
        "--recursive",
        output_directory,
        os.path.join(s3_path, output_name)
    ]

    if dry:
        aws_args = ["echo"] + aws_args

    run(aws_args, check=True)


def main(shasta_path, fastq_paths, s3_path, dry):
    for path in fastq_paths.split(','):
        run_evaluation(shasta_path=shasta_path, fastq_path=path, s3_path=s3_path, dry=dry)
        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--shasta","-s",
        type=str,
        required=True,
        help="Path of shasta binary file to use. Can also be a comma separated list of paths"
    )

    parser.add_argument(
        "--input","-i",
        type=str,
        required=True,
        help="path of fastq file to be evaluated"
    )

    parser.add_argument(
        "--s3",
        type=str,
        required=True,
        help="AWS s3 URI of directory to upload to"
    )

    parser.add_argument(
        "--dry",
        dest="dry",
        action="store_true",
        required=False,
        help="echo the commands instead of executing them"
    )

    args = parser.parse_args()

    main(
        shasta_path=args.shasta,
        fastq_paths=args.input,
        s3_path=args.s3,
        dry=args.dry,
    )
