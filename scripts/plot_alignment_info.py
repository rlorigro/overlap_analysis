from matplotlib import pyplot
import numpy
import sys


def load_dataset_as_numpy(csv_path):
    """
    csv has columns: readId1, readId2, minAlignedFraction, markerCount, maxDrift, maxSkip, trim, overlap
    :param csv_path:
    :return:
    """
    data = list()

    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            tokens = line.strip().split(',')
            line_data = list(map(float, tokens))

            data.append(line_data)

    data = numpy.array(data)

    return data


def main(csv_path):
    data = load_dataset_as_numpy(csv_path=csv_path)

    # data = data[:10000,:]

    fig, axes = pyplot.subplots(nrows=5,ncols=5)

    label = data[:,-1]
    label = (label > 0)
    inverse_label = numpy.invert(label)

    true_color = [0.129,0.216,0.667]
    false_color = [0.004,0.569,0.216]

    labels = ["minAlignedFraction", "markerCount", "maxDrift", "maxSkip", "trim"]

    for i in range(5):
        for j in range(5):
            if (i) + (4 - j) < 4:
                axes[i][j].axis("off")
                continue

            if i == j:
                y = data[:,2+i]
                y_false = y[inverse_label]
                y_true = y[label]

                start = numpy.min(y)
                stop = numpy.max(y)

                n_bins = None
                if i > 1:
                    n_bins = min(stop - start, 500)
                else:
                    n_bins = 500

                step = float(stop - start)/n_bins

                bins = numpy.arange(start, stop + step, step=step)
                y_true_frequencies, _ = numpy.histogram(y_true, bins=bins)
                y_false_frequencies, _ = numpy.histogram(y_false, bins=bins)

                # y_true_frequencies = numpy.log10(y_true_frequencies)
                # y_false_frequencies = numpy.log10(y_false_frequencies)

                x = bins[:-1] + (step/2)

                axes[i][j].plot(x, y_true_frequencies, color=true_color, linewidth=0.2)
                axes[i][j].plot(x, y_false_frequencies, color=false_color, linewidth=0.2)
                axes[i][j].fill_between(x=x, y1=y_true_frequencies, y2=0, color=true_color, alpha=0.3)
                axes[i][j].fill_between(x=x, y1=y_false_frequencies, y2=0, color=false_color, alpha=0.3)

                axes[i][0].set_ylabel(labels[i])
                axes[-1][j].set_xlabel(labels[j])

                if labels[i] == labels[j] == "markerCount":
                    axes[i][j].set_xlim([-200,10000])

            else:
                x = data[:,2+j]
                y = data[:,2+i]

                x_false = x[inverse_label]
                y_false = y[inverse_label]
                x_true = x[label]
                y_true = y[label]

                axes[i][j].scatter(x=x_true,y=y_true, color=true_color, s=0.2, alpha=0.08)
                axes[i][j].scatter(x=x_false,y=y_false, color=false_color, s=0.2, alpha=0.08)

                if labels[i] == "markerCount":
                    axes[i][j].set_ylim([-200,10000])

                elif labels[j] == "markerCount":
                    axes[i][j].set_xlim([-200,10000])

    pyplot.show()
    pyplot.close()


if __name__ == "__main__":

    if len(sys.argv) != 2:
        exit("ERROR: need to provide 1 argument: path of input GFA")

    main(sys.argv[1])
