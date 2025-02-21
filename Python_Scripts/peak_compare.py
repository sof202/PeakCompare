import argparse
from create_confidence_intervals import generate_pvalue_ci
from determine_metric import calculate_metric
from determine_psuedo_peaks import compare_pvalue_ci, determine_psuedopeaks
from extract_region import extract_bedbase_region
from label_peak_type import label_peak_type, convert_narrow_peak_to_bedbase
from IO import BedGraph, Bed
from multiprocessing import Pool
import pandas as pd
import sys


def string_list(arg):
    return [x for x in arg.split(',')]


def integer_list(arg):
    try:
        return [int(x) for x in arg.split(',')]
    except ValueError:
        raise argparse.ArgumentTypeError(
            "Argument must be a comma-separated list of integers")


def run(reference_merged_peaks: Bed,
        reference_unmerged_peaks: Bed,
        reference_bias_track: BedGraph,
        reference_coverage_track: BedGraph,
        comparison_bias_track: BedGraph,
        comparison_coverage_track: BedGraph,
        comparison_pvalue_track: BedGraph,
        significance: float,
        window_size: int,
        cutoff: float,
        unmerged: bool,
        chromosome: str,
        start: int,
        end: int) -> float:
    reference_merged_peaks = convert_narrow_peak_to_bedbase(
        reference_merged_peaks,
        chromosome,
        start,
        end
    )
    reference_unmerged_peaks = convert_narrow_peak_to_bedbase(
        reference_unmerged_peaks,
        chromosome,
        start,
        end
    )
    reference_bias_track = extract_bedbase_region(
        reference_bias_track,
        chromosome,
        start,
        end
    )
    reference_coverage_track = extract_bedbase_region(
        reference_coverage_track,
        chromosome,
        start,
        end
    )
    comparison_bias_track = extract_bedbase_region(
        comparison_bias_track,
        chromosome,
        start,
        end
    )
    comparison_coverage_track = extract_bedbase_region(
        comparison_coverage_track,
        chromosome,
        start,
        end
    )
    comparison_pvalue_track = extract_bedbase_region(
        comparison_pvalue_track,
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
        significance,
        window_size
    )
    comparison_pvalue_ci = generate_pvalue_ci(
        comparison_bias_track,
        comparison_coverage_track,
        significance,
        window_size
    )
    compared_pvalues = compare_pvalue_ci(
        reference_pvalue_ci,
        comparison_pvalue_ci
    )
    pseudopeaks = determine_psuedopeaks(
        comparison_pvalue_track,
        compared_pvalues,
        reference_labelled_peaks,
        cutoff
    )
    metric = calculate_metric(
        reference_labelled_peaks,
        pseudopeaks,
        include_merged_peaks=(not unmerged)
    )
    return chromosome, start, end, metric


def main(args: argparse.Namespace, regions: list) -> None:
    reference_merged_peaks = Bed.read_from_file(
        args.reference_merged_peaks_file)
    reference_unmerged_peaks = Bed.read_from_file(
        args.reference_unmerged_peaks_file)
    reference_bias_track = BedGraph.read_from_file(
        args.reference_bias_track_file)
    reference_coverage_track = BedGraph.read_from_file(
        args.reference_coverage_track_file)
    comparison_bias_track = BedGraph.read_from_file(
        args.comparison_bias_track_file)
    comparison_coverage_track = BedGraph.read_from_file(
        args.comparison_coverage_track_file)
    comparison_pvalue_track = BedGraph.read_from_file(
        args.comparison_pvalue_file)
    run_arguments = [
        (
            reference_merged_peaks, reference_unmerged_peaks,
            reference_bias_track, reference_coverage_track,
            comparison_bias_track, comparison_coverage_track,
            comparison_pvalue_track, args.significance, args.window_size,
            args.cutoff, args.unmerged, chromosome, start, end
        )
        for chromosome, start, end in zip(
            args.chromosome, args.start, args.end
        )
    ]
    with Pool() as pool:
        results = pd.DataFrame(pool.starmap(run, run_arguments))
    results.to_csv(sys.stdout, header=False, index=False,
                   float_format='%.4f', sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="PeakCompare",
        description="Determine how likely a peak is to be in two datasets."
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
        type=string_list,
        help=("The chromosomes of the regions you wish to inspect. Comma "
              "separated list.")
    )
    parser.add_argument(
        "start",
        type=integer_list,
        help=("The base pair positions at the start of the regions you wish "
              "to inspect. Comma separated list.")
    )
    parser.add_argument(
        "end",
        type=integer_list,
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
    main(args)
