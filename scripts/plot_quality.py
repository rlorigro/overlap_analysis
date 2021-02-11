from matplotlib import pyplot
import ruptures
import scipy.stats
import random
import numpy
import math
import sys


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

                # print(qualities_string[:10])
                # print(q_values[:10])
                # print(qualities[:10])
                # print(smoothed_qualities[:10])

                print(n, column)

                # max_breakpoint = optimal_partition(smoothed_qualities)
                # print("Breakpoint: ", max_breakpoint)

                model = "l1"  # "l2", "rbf"

                # PELT
                # algo = ruptures.Pelt(model=model, min_size=100, jump=250).fit(smoothed_qualities)
                # breakpoints = algo.predict(pen=200)

                # BinSeg
                # dimensions = 1      # Dimensionality of data
                # sigma = 0.01        # Noise level (stddev)
                # penalty = numpy.log(n)*dimensions*sigma**2
                # algo = ruptures.Binseg(model=model, min_size=1000, jump=250).fit(smoothed_qualities)
                # breakpoints = algo.predict(pen=200)

                # print(breakpoints)

                axis = None
                if column is None:
                    axis = axes[n-1]
                else:
                    axis = axes[n-1][column]

                axis.plot(smoothed_qualities)

                # axis.legend([name])

                axis.set_ylim([0,0.5])

                # if max_breakpoint.t is not None:
                #     axis.axvline(max_breakpoint.t, color="orange")

                # for x in breakpoints:
                #     axis.axvline(x, color="orange")

                if n == n_reads:
                    break


def log_likelihood(x, t1, t2):

    mean = numpy.mean(x[t1:t2])

    # print("mean (%d,%d):" % (t1,t2), mean)

    ll_a = numpy.log(2*numpy.pi)
    ll_b = (x[t1:t2]-mean)**2

    ll = (t2-t1)*(ll_a + numpy.log(numpy.sum(ll_b)/(t2-t1) + 1))

    return ll


class Breakpoint:
    def __init__(self, t, score):
        self.t = t
        self.score = score

    def __str__(self):
        return str(self.t) + " " + str(self.score)


def scan_breakpoints(x, step_size, beta):
    results = list()

    for t in range(step_size, len(x), step_size):
        ll_left = log_likelihood(x, t1=0, t2=t)
        ll_right = log_likelihood(x, t1=t, t2=len(x))
        ll_all = log_likelihood(x, t1=0, t2=len(x))

        # print(ll_left, ll_right)
        # print(ll_left + ll_right, ll_all)
        # print(ll_all - (ll_left + ll_right))
        # print()

        score = ll_all - (ll_left + ll_right) - beta

        b = Breakpoint(t=t, score=score)

        results.append(b)

    return results


def optimal_partition(x):
    beta = 10
    step_size = 250

    breakpoints = scan_breakpoints(x, beta=beta, step_size=step_size)

    max_breakpoint = Breakpoint(t=None, score=-sys.maxsize)

    for item in breakpoints:
        if item.score > 0 and item.score > max_breakpoint.score:
            max_breakpoint = item

    return max_breakpoint


def segment_by_t_test(x):
    midpoint = int(len(x)/2)

    n_left = midpoint
    n_right = len(x) - midpoint

    samples_left = min(1000, n_left)
    samples_right = min(1000, n_right)

    x_left = numpy.random.choice(x[:midpoint], size=samples_left)
    x_right = numpy.random.choice(x[midpoint:], size=samples_right)

    mean_left = numpy.mean(x_left)
    mean_right = numpy.mean(x_right)

    variance_left = numpy.var(x_left)
    variance_right = numpy.var(x_right)

    print("n:\t\t", n_left, n_right)
    print("samples:\t", samples_left, samples_right)
    print("mean:\t\t", mean_left, mean_right)
    print("variance:\t", variance_left, variance_right)

    t = (mean_left - mean_right) / math.sqrt(variance_left/samples_left + variance_right/samples_right)

    print("t:\t\t",t)

    v_numerator = ((variance_left/samples_left) + (variance_right/samples_right))**2
    v_denominator = ((variance_left/samples_left)**2)/(samples_left-1)+((variance_right/samples_right)**2)/(samples_right-1)

    v = v_numerator / v_denominator

    print("v:\t\t",v)

    # p_a = math.gamma((v+1)/2)/(math.sqrt(v*math.pi)*math.gamma(v/2))
    # p_b = (1+(t**2)/v)**((v**2)/2)
    #
    # p = p_a*p_b

    t_dist = scipy.stats.t(v)
    p = t_dist.cdf(t)

    print("p:\t\t",p)

    a, p2 = scipy.stats.ttest_ind(x_left, x_right, equal_var=False)

    print("p2:\t\t", p2)


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

    n_reads = max(n_reads_a, n_reads_b)

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


if __name__ == "__main__":
    # test_convolution()
    main()
