import argparse
from create_confidence_intervals import generate_pvalue_ci
from determine_metric import calculate_metric
from determine_psuedo_peaks import compare_pvalue_ci, determine_psuedopeaks
from extract_region import extract_bedbase_region
from label_peak_type import label_peak_type, convert_narrow_peak_to_bedbase
from IO import BedGraph, Bed


def list_type(arg):
    try:
        return [int(x) for x in arg.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Argument must be a comma-separated list of integers")


def process_regions(args: argparse.Namespace) -> list:
    regions = [(chromosome, start, end) for chromosome, start,
               end in zip(args.chromosome, args.start, args.end)]
    return regions


def main(args: argparse.Namespace) -> None:
    chromosome = args.chromosome
    start = args.start
    end = args.end

    reference_merged_peaks = convert_narrow_peak_to_bedbase(
        Bed.read_from_file(args.reference_merged_peaks_file),
        chromosome,
        start,
        end
    )
    reference_unmerged_peaks = convert_narrow_peak_to_bedbase(
        Bed.read_from_file(args.reference_unmerged_peaks_file),
        chromosome,
        start,
        end
    )
    reference_bias_track = extract_bedbase_region(
        BedGraph.read_from_file(args.reference_bias_track_file),
        chromosome,
        start,
        end
    )
    reference_coverage_track = extract_bedbase_region(
        BedGraph.read_from_file(args.reference_coverage_track_file),
        chromosome,
        start,
        end
    )
    comparison_bias_track = extract_bedbase_region(
        BedGraph.read_from_file(args.comparison_bias_track_file),
        chromosome,
        start,
        end
    )
    comparison_coverage_track = extract_bedbase_region(
        BedGraph.read_from_file(args.comparison_coverage_track_file),
        chromosome,
        start,
        end
    )
    comparison_pvalue_track = extract_bedbase_region(
        BedGraph.read_from_file(args.comparison_pvalue_file),
        chromosome,
        start,
        end
    )
    reference_labelled_peaks = label_peak_type(
        reference_merged_peaks,
        reference_unmerged_peaks
    )
    reference_pvalue_ci = generate_pvalue_ci(
        reference_bias_track,
        reference_coverage_track,
        args.significance,
        args.window_size
    )
    comparison_pvalue_ci = generate_pvalue_ci(
        comparison_bias_track,
        comparison_coverage_track,
        args.significance,
        args.window_size
    )
    compared_pvalues = compare_pvalue_ci(
        reference_pvalue_ci,
        comparison_pvalue_ci
    )
    pseudopeaks = determine_psuedopeaks(
        comparison_pvalue_track,
        compared_pvalues,
        reference_labelled_peaks,
        args.cutoff
    )
    metric = calculate_metric(
        reference_labelled_peaks,
        pseudopeaks,
        include_merged_peaks=(not args.unmerged)
    )
    if args.parsable:
        print(metric)
    else:
        print("The metric for these datasets over the range ",
              f"{start} to {end} for chromosome {chromosome} is: {metric}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="PeakCompare",
        description="Determine how likely a peak is to be in two datasets."
    )
    parser.add_argument(
        "-p",
        "--parsable",
        action="store_true",
        help=("Set this if you want the output of the script to be computer ",
              "parsable (and not human readable).")
    )
    parser.add_argument(
        "--unmerged",
        action="store_true",
        help=("Set this if you want to discount peaks that are a result of "
              "merging when calculating the metric.")
    )
    parser.add_argument(
        "--significance",
        nargs='?',
        const=0.95,
        default=0.95,
        type=float,
        help=("The significance used when calculating confidence intervals.")
    )
    parser.add_argument(
        "--window_size",
        nargs='?',
        const=50,
        default=50,
        type=int,
        help=("The window size used when calculating confidence intervals.")
    )
    parser.add_argument(
        "chromosome",
        type=list_type,
        help=("The chromosomes of the regions you wish to inspect. Comma "
              "separated list.")
    )
    parser.add_argument(
        "start",
        type=list_type,
        help=("The base pair positions at the start of the regions you wish "
              "to inspect. Comma separated list.")
    )
    parser.add_argument(
        "end",
        type=list_type,
        help=("The base pair positions at the end of the regions you wish to "
              "inspect. Comma separated list")
    )
    parser.add_argument(
        "reference_merged_peaks_file",
        help=("The narrow peak file from reference dataset where peaks are"
              "merged.")
    )
    parser.add_argument(
        "reference_unmerged_peaks_file",
        help=("The narrow peak file from reference dataset where peaks are"
              "not merged.")
    )
    parser.add_argument(
        "reference_bias_track_file",
        help="The bias track file for the reference dataset."
    )
    parser.add_argument(
        "reference_coverage_track_file",
        help="The coverage track (pileup) file for the reference dataset."
    )
    parser.add_argument(
        "comparison_bias_track_file",
        help="The bias track file for the comparison dataset."
    )
    parser.add_argument(
        "comparison_coverage_track_file",
        help="The coverage track (pileup) file for the comparison dataset."
    )
    parser.add_argument(
        "comparison_pvalue_file",
        help="The pvalues for the comparison dataset."
    )
    parser.add_argument(
        "cutoff",
        type=float,
        help="The cutoff used to call peaks in the reference dataset."
    )
    args = parser.parse_args()
    regions = process_regions(args)
    main(args)
