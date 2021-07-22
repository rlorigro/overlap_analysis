from modules.Fastx import build_index


def main():
    fastq_path = "/home/ryan/data/nanopore/human/test/cross_strand_analysis/whole_flowcell/palindromes.fastq"
    faidx_path = build_index(fastq_path)

    thresholds = [0,5_000,10_000,50_000,100_000]
    total_lengths = [0 for i in thresholds]

    with open(faidx_path, 'r') as file:
        for l,line in enumerate(file):
            data = line.strip().split()

            length = int(data[1])

            for i,t in enumerate(thresholds):
                if length > t:
                    total_lengths[i] += length

    for i,t in enumerate(thresholds):
        print("length greater than %d:\t%d" % (t,total_lengths[i]))





if __name__ == "__main__":
    main()
