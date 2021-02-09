from matplotlib import pyplot
import numpy


def plot_fastq(fastq_path, axes, n_reads, column=None):
    window_size = 200
    step_size = 6.0/float(window_size)

    x_kernel = numpy.arange(start=-3,stop=3+step_size,step=step_size)
    kernel = (2.71828**(-(x_kernel**2)/2))/numpy.sqrt(2*3.14159)
    kernel /= numpy.sum(kernel)

    n = 0
    with open(fastq_path, 'r') as file:
        name = None
        qualities = None

        for l,line in enumerate(file):
            if l % 4 == 0:
                n += 1
                name = line.strip()
            elif l % 4 == 3:
                qualities_string = line.strip()

                q_values = [(q-33) for q in map(ord, list(qualities_string))]
                qualities = numpy.array([(10**(float(q-33)/-10)) for q in map(ord, list(qualities_string))], dtype=numpy.float)

                smoothed_qualities = numpy.convolve(qualities, kernel, 'valid')

                print(qualities_string[:10])
                print(q_values[:10])
                print(qualities[:10])
                print(smoothed_qualities[:10])

                axis = None
                if column is None:
                    axis = axes[n-1]
                else:
                    axis = axes[n-1][column]

                axis.plot(smoothed_qualities)

                axis.set_ylim([0,0.5])

                if n == n_reads:
                    break


def main():
    cross_strand_fastq_path = "/home/ryan/data/nanopore/human/test/cross_strand_analysis/GM24385_2_Guppy_3.6.0_prom_cross_strand_reads.fastq"
    random_fastq_path = "/home/ryan/data/nanopore/human/test/cross_strand_analysis/GM24385_2_Guppy_3.6.0_prom_random_reads.fastq"

    n_reads_a = 0
    with open(cross_strand_fastq_path, 'r') as file:
        for l,line in enumerate(file):
            if l % 4 == 0:
                n_reads_a += 1

    n_reads_b = 0
    with open(random_fastq_path, 'r') as file:
        for l,line in enumerate(file):
            if l % 4 == 0:
                n_reads_b += 1

    n_reads = min(n_reads_a, n_reads_b)

    figure,axes = pyplot.subplots(nrows=n_reads, ncols=2)
    figure.set_size_inches(12,6)
    pyplot.subplots_adjust(hspace=.7)

    plot_fastq(fastq_path=cross_strand_fastq_path, axes=axes, n_reads=n_reads, column=0)
    plot_fastq(fastq_path=random_fastq_path, axes=axes, n_reads=n_reads, column=1)

    axes[-1][0].set_xlabel("Sliding Window (200bp Gaussian)")
    axes[-1][1].set_xlabel("Sliding Window (200bp Gaussian)")

    axes[int(n_reads/2)][0].set_ylabel("P(error)")
    axes[int(n_reads/2)][1].set_ylabel("P(error)")

    axes[0][0].set_title("Cross-strand Reads")
    axes[0][1].set_title("Random Reads")

    pyplot.show()
    pyplot.close()


def segment(x):
    midpoint = len(x)/2






if __name__ == "__main__":
    # test_convolution()
    main()
