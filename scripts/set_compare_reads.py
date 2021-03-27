import argparse
import sys
import os


def parse_read_id_from_file(path):
    ids = set()

    with open(path, 'r') as file:
        for l,line in enumerate(file):
            line = line.strip()
            ids.add(line)

    return ids


def main(a_path, b_path):
    a = parse_read_id_from_file(a_path)
    b = parse_read_id_from_file(b_path)

    a_only = a - b
    b_only = b - a
    a_and_b = a.intersection(b)

    print("A only: %d" % len(a_only))
    # for item in a_only:
    #     print(item)

    print("B only: %d" % len(b_only))
    # for item in b_only:
    #     print(item)

    print("A and B: %d" % len(a_and_b))
    # for item in a_and_b:
    #     print(item)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-a",
        type=str,
        required=True,
        help="path of file containing a list of read ids, one on each line"
    )
    parser.add_argument(
        "-b",
        type=str,
        required=True,
        help="path of file containing a list of read ids, one on each line"
    )

    args = parser.parse_args()

    main(
        a_path=args.a,
        b_path=args.b
    )
