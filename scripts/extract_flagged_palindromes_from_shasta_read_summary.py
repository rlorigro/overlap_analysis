

def extract_read_ids_using_read_names_from_shasta_read_summary(summary_path, read_names):
    read_names = {r for r in read_names}

    with open(summary_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split(',')

            id = data[0]
            name = data[1]

            if name in read_names:
                print(id + ',' + name)
                read_names.remove(name)

    return


"""
Id,Name,RawLength,RleLength,RawOverRleLengthRatio,MarkerCount,MarkerDensity,MaximumMarkerOffset,Palindromic,Chimeric,AlignmentCandidates,ReadGraphNeighbors,VertexCount,VertexDensity,runid,sampleid,read,ch,start_time,
"""
def main():
    summary_path = "/home/ryan/data/nanopore/human/test/overlap_guppy_360_whole_flowcell/run3/ReadSummary.csv"
    other_summary_path = "/home/ryan/data/nanopore/human/test/overlap_guppy_360_whole_flowcell/run1/ReadSummary.csv"

    n_palindromic = 0

    palindromic_names = list()

    with open(summary_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split(',')

            id = data[0]
            name = data[1]

            palindromic = False
            if data[8] == "Yes":
                palindromic = True
            elif data[8] == "No":
                palindromic = False
            else:
                exit("ERROR: unrecognized value for palindromic flag on line " + str(l) + " : " + line)

            if palindromic:
                # print(name, id)
                n_palindromic += 1
                palindromic_names.append(name)

    # print("n_palindromic: ", n_palindromic)

    extract_read_ids_using_read_names_from_shasta_read_summary(other_summary_path, read_names=palindromic_names)


if __name__ == "__main__":
    main()

