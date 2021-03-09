from matplotlib import pyplot
from modules.IterativeHistogram import IterativeHistogram
from modules.Fastq import *
import argparse
import mmap
import os


def plot_length_distribution(histogram, axes, label, weighted=False, normalize=False):

    if normalize:
        frequencies = histogram.get_normalized_histogram()
    else:
        frequencies = histogram.get_histogram()

    bin_centers = histogram.get_bin_centers()

    if weighted:
        frequencies *= bin_centers

    axes.plot(bin_centers, frequencies, alpha=0.7, label=label)

    axes.set_xlabel("Length (bp)")

    if weighted:
        axes.set_ylabel("Bases")
    elif not normalize:
        axes.set_ylabel("Frequency")
    else:
        axes.set_ylabel("Density")


def compute_subsample_rates(template_distribution, subsample_distribution, n_bins, max_reduction_rate):
    bin_edges = template_distribution.get_bin_edges()

    subsample_rates = [0]*(n_bins)

    sys.stderr.write("i\tbin\tratio\n")

    normalized_template_distribution = template_distribution.get_normalized_histogram()
    normalized_subsample_distribution = subsample_distribution.get_normalized_histogram()

    for i in range(n_bins):
        ratio = 0
        if normalized_template_distribution[i] > 0 and normalized_subsample_distribution[i] > 0:
            ratio = float(normalized_template_distribution[i])/normalized_subsample_distribution[i]

        subsample_rates[i] = ratio

        sys.stderr.write("%d\t%.2f\t%.2f\n" % (i, bin_edges[i], ratio))

    bin_indexes = [i for i in range(len(subsample_rates))]
    subsample_coverages = subsample_distribution.get_normalized_histogram() * subsample_distribution.get_bin_centers()
    template_coverages = template_distribution.get_normalized_histogram() * template_distribution.get_bin_centers()

    worst_clipped_indexes = sorted(bin_indexes, key=lambda x: template_coverages[x] - subsample_coverages[x], reverse=True)

    # Compute the reduction ratio from an average of the worst bins with a ratio >1 (most clipping)
    reduction_rate = [i for i in worst_clipped_indexes[:20] if i > 1]
    reduction_rate = sum(reduction_rate)/len(reduction_rate) if len(reduction_rate) > 0 else 1

    reduction_rate = min(reduction_rate, max_reduction_rate)

    # Scale down the ratio so that even if there are subsample:template ratios >1, the overall distribution can match
    # by shrinking the total coverage
    subsample_rates = [r/reduction_rate for r in subsample_rates]

    return subsample_rates


def main(template_path, subsample_path, output_dir, max_bases, min_length, max_reduction_rate, dry_run):
    if not os.path.exists(output_dir):
        sys.stderr.write("Creating output dir: %s\n" % output_dir)
        os.makedirs(output_dir)

    output_path = os.path.join(output_dir, "output.fastq")

    # Allow an fai to be used instead of a full fastq for the template
    if not template_path.endswith(".fai"):
        template_faidx_path = build_index(template_path)
    else:
        template_faidx_path = template_path

    subsample_faidx_path = build_index(subsample_path)

    template_name_to_offset, template_index_elements = load_fastq_index(faidx_path=template_faidx_path)
    subsample_name_to_offset, subsample_index_elements = load_fastq_index(faidx_path=subsample_faidx_path)

    n_bins = 500
    start = 0
    stop = 500_000

    template_distribution = IterativeHistogram(start=start, stop=stop, n_bins=n_bins)
    subsample_distribution = IterativeHistogram(start=start, stop=stop, n_bins=n_bins)
    resulting_distribution = IterativeHistogram(start=start, stop=stop, n_bins=n_bins)

    sys.stderr.write("Counting template and subsample distributions...\n")
    sys.stderr.flush()

    # Find the length distribution for the template file
    for item in template_index_elements:
        if item.length > min_length:
            template_distribution.update(item.length)

    # Find the length distribution for the file to be subsampled
    for item in subsample_index_elements:
        if item.length > min_length:
            subsample_distribution.update(item.length)

    # Default to the full coverage of the subsample file (as the max coverage, often not achieved)
    if max_bases is None:
        max_bases = sum(subsample_distribution.get_histogram()*subsample_distribution.get_bin_centers())

    subsample_rates = compute_subsample_rates(template_distribution, subsample_distribution, n_bins, max_reduction_rate)

    sys.stderr.write("Subsampling FASTQ...\n")
    subsampled_bases = 0
    with open(subsample_path, 'rb') as input_file, open(output_path, 'wb') as output_file:
        mm = mmap.mmap(input_file.fileno(), 0, prot=mmap.PROT_READ)

        for name, i in subsample_name_to_offset.items():
            if i % 1000 == 0:
                sys.stderr.write("\r{:2.1%}".format(i / len(subsample_name_to_offset)))

            offset_index = subsample_name_to_offset[name]
            index_element = subsample_index_elements[offset_index]

            # Skip reads below a threshold (nanopore has many reads < 100bp)
            if index_element.length < min_length:
                continue

            bin_index = subsample_distribution.get_bin(index_element.length)

            # Ignore reads longer than the histogram bounds
            if bin_index >= n_bins:
                continue

            cutoff = subsample_rates[bin_index] * 101

            # Perform a simple mod-based random sampling to decide whether to pick this read
            if i % 101 > cutoff:
                continue

            subsampled_bases += index_element.length

            resulting_distribution.update(index_element.length)

            if not dry_run:
                s = extract_bytes_from_file(mmap_file_object=mm,
                                            offset=index_element.sequence_offset,
                                            n_bytes=index_element.length)

                q = extract_bytes_from_file(mmap_file_object=mm,
                                            offset=index_element.quality_offset,
                                            n_bytes=index_element.length)

                output_file.write(b'@')
                output_file.write(name.encode('utf-8'))
                output_file.write(b'\n')
                output_file.write(s)
                output_file.write(b'\n')
                output_file.write(b'+')
                output_file.write(b'\n')
                output_file.write(q)
                output_file.write(b'\n')

            if subsampled_bases >= max_bases:
                sys.stderr.write("Stopping because bases sampled (%d) >= max_bases (%d)\n" % (subsampled_bases, max_bases))
                break

    sys.stderr.write('\n')

    figure = pyplot.figure()
    axes = pyplot.axes()

    image_path = os.path.join(output_dir, "normalized_distributions.png")
    plot_length_distribution(histogram=template_distribution, axes=axes, label="template", normalize=True)
    plot_length_distribution(histogram=subsample_distribution, axes=axes, label="subsample", normalize=True)
    plot_length_distribution(histogram=resulting_distribution, axes=axes, label="resulting", normalize=True)

    figure.set_size_inches(16,8)
    axes.legend()
    sys.stderr.write("Saving plot to: %s\n" % image_path)
    pyplot.savefig(image_path, dpi=200)

    figure2 = pyplot.figure()
    axes2 = pyplot.axes()

    image_path = os.path.join(output_dir, "coverage_weighted_distributions.png")
    plot_length_distribution(histogram=template_distribution, weighted=True, axes=axes2, label="template")
    plot_length_distribution(histogram=subsample_distribution, weighted=True, axes=axes2, label="subsample")
    plot_length_distribution(histogram=resulting_distribution, weighted=True, axes=axes2, label="resulting")

    figure2.set_size_inches(16,8)
    axes2.legend()
    sys.stderr.write("Saving plot to: %s\n" % image_path)
    pyplot.savefig(image_path, dpi=200)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--template",
        type=str,
        required=True,
        help="path of FASTQ file, OR .fai file from a FASTQ (attempt to match this distribution)"
    )
    parser.add_argument(
        "--subsample",
        type=str,
        required=True,
        help="path of file containing FASTQ sequence"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="path to directory to dump output"
    )
    parser.add_argument(
        "--max_bases",
        type=int,
        required=False,
        default=None,
        help="maximum amount of bases to subsample"
    )
    parser.add_argument(
        "--min_length",
        type=int,
        required=False,
        default=100,
        help="ignore all reads below this length"
    )
    parser.add_argument(
        "--max_reduction_rate","-r",
        type=float,
        required=False,
        default=1.5,
        help="To accommodate certain distributions, it may require scaling down the entire dataset. Set this value to 1 "
             "to completely prevent downscaling, or to choose a number (1,inf) to set the upper limit on reduction rate. "
             "Example: r=1.5, then the total coverage may be allowed to reduce by 1/1.5 = 1/3"
    )
    parser.add_argument(
        "--dry_run","-d",
        dest="dry_run",
        required=False,
        action="store_true",
        help="Add this boolean flag to compute input and result distributions without actually reading or writing any "
             "fastq (only reads the associated .fai index files)"
    )

    args = parser.parse_args()

    main(
        template_path=args.template,
        subsample_path=args.subsample,
        output_dir=args.output_dir,
        max_bases=args.max_bases,
        min_length=args.min_length,
        max_reduction_rate=args.max_reduction_rate,
        dry_run=args.dry_run
    )


