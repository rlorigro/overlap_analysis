from collections import defaultdict
from sklearn import metrics
import argparse
import os


class PhasedRead:
    def __init__(self, id, name, component, phase):
        self.id = id
        self.name = name
        self.component = int(component)
        self.phase = int(phase)


# Iterate the shasta CSV and yield one PhasedRead data object at a time
def parse_shasta_phase_csv(shasta_path):
    with open(shasta_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                # Skip header line
                continue

            # OrientedReadId,ReadName,Component,Phase
            tokens = line.strip().split(',')

            # Ignore unphased reads (for now)
            if (len(tokens[-1]) == 0) and (len(tokens[-2]) == 0):
                continue

            yield PhasedRead(tokens[0], tokens[1], tokens[2], tokens[3])


def main(hap1_margin_path, hap2_margin_path, shasta_path):
    discordant_csv_path = os.path.join(os.path.dirname(shasta_path), "discordantPhasedReads.csv")

    margin_phases = (set(),set())
    shasta_components = defaultdict(lambda: (set(),set()))
    name_to_id = dict()

    # Read all the Shasta phase assignments
    for read in parse_shasta_phase_csv(shasta_path):
        shasta_components[read.component][read.phase].add(read.name)
        name_to_id[read.name] = read.id

        if read.name in shasta_components[read.component][1-read.phase]:
            exit("ERROR: read contained both in hap1 and hap2 of Shasta phases: " + read.name)

    # Read the Margin phase assignments
    with open(hap1_margin_path, 'r') as margin0_file, open(hap2_margin_path, 'r') as margin1_file:
        for line in margin0_file:
            name = line.strip()
            margin_phases[0].add(name)

        for line in margin1_file:
            name = line.strip()
            margin_phases[1].add(name)

    # Compare SHasta to Margin and write results to CSV + stdout
    with open(discordant_csv_path, 'w') as file:
        file.write("OrientedReadId,ReadName,Component,ShastaPhase,MarginPhase")

        for component_id,component in shasta_components.items():
            print("Evaluating component:", component_id)

            names = list()
            shasta_phase_array = list()
            margin_phase_array = list()

            # Compare phase assignments for Shasta and the user-provided dataset, for each read in each component
            for p,phase in enumerate(component):
                for read_name in phase:
                    margin_phase = -1

                    if read_name in margin_phases[0]:
                        margin_phase = 0

                    if read_name in margin_phases[1]:
                        if margin_phase == 0:
                            exit("ERROR: read contained both in hap1 and hap2 of Margin phases: " + read_name)
                        else:
                            margin_phase = 1

                    # Only store results for reads that also exist in the user-provided dataset
                    if margin_phase >= 0:
                        names.append(read_name)
                        shasta_phase_array.append(p)
                        margin_phase_array.append(margin_phase)

            count_matrix = [[0,0],[0,0]]
            confusion_matrix = [[0,0],[0,0]]

            for i in range(len(names)):
                count_matrix[0][shasta_phase_array[i]] += 1
                count_matrix[1][margin_phase_array[i]] += 1

                confusion_matrix[shasta_phase_array[i]][margin_phase_array[i]] += 1

            # Compute the weight of the diagonal and antidiagonal, and infer which is True Positives
            diagonal_sum = confusion_matrix[0][0] + confusion_matrix[1][1] + 1e-12
            antidiagonal_sum = confusion_matrix[1][0] + confusion_matrix[0][1]

            if 2 > (antidiagonal_sum / diagonal_sum) > 0.5:
                print("WARNING: ratio of antidiagonal_sum:diagonal_sum is only " + str(antidiagonal_sum / diagonal_sum))

            # Indicate which coordinates in the confusion matrix represent discordant reads between shasta's
            # and the user-provided phase assignments
            discoordinate_a = (0,1)
            discoordinate_b = (1,0)

            # Iterate reads again and write the ones that are discordant with the user-provided phases
            if antidiagonal_sum > diagonal_sum:
                discoordinate_a = (0,0)
                discoordinate_b = (1,1)

            for i in range(len(names)):
                coordinate = (shasta_phase_array[i], margin_phase_array[i])

                if coordinate == discoordinate_a or coordinate == discoordinate_b:
                    file.write(name_to_id[names[i]])
                    file.write(',')
                    file.write(names[i])
                    file.write(',')
                    file.write(str(component_id))
                    file.write(',')
                    file.write(str(shasta_phase_array[i]))
                    file.write(',')
                    file.write(str(margin_phase_array[i]))
                    file.write('\n')

            print("Shasta phase counts:",count_matrix[0])
            print("Margin phase counts:",count_matrix[1])

            print("Confusion matrix:")
            print(confusion_matrix[0])
            print(confusion_matrix[1])

            adjusted_rand_index = metrics.adjusted_rand_score(margin_phase_array, shasta_phase_array)
            print("ARI:", adjusted_rand_index)
            print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-1","--hap1",
        type=str,
        required=True,
        help="path of marginphase file for hap1"
    )
    parser.add_argument(
        "-2","--hap2",
        type=str,
        required=True,
        help="path of marginphase file for hap2"
    )
    parser.add_argument(
        "-s","--shasta",
        type=str,
        required=True,
        help="path of shasta phase components csv"
    )

    args = parser.parse_args()

    main(
        hap1_margin_path=args.hap1,
        hap2_margin_path=args.hap2,
        shasta_path=args.shasta,
    )

