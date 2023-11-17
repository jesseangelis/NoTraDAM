import argparse
from src.main import Pipeline

def get_arguments():
    """
    Parses the arguments

    Returns
    -----------
    args : argparse.ArgumentParser object
        arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference_annotation", type=str, required=True, help="The reference annotation as .gtf file")
    parser.add_argument("-n", "--query_annotations", nargs="+", required=True, help="The predicted annotation or annotations as .gtf file/s")
    parser.add_argument("-i", "--transcript_ids", type=str, required=False, help="Text file, containing the transcript IDs, that were removed from the reference annotation. If not given, all transcripts are considred removed")
    parser.add_argument("-f", "--genome_fasta", type=str, required=True, help="Fasta file with chromosome sequences. The chromosome IDs must be the same as in the annotations")
    parser.add_argument("-o", "--output_directory", type=str, required=False, help="Directory used for the outputs")
    parser.add_argument("-k", "--keep_temporary", action="store_true", required=False, help="If given, temporay files are kept. This may increase perfrmances for reruns, but may result in wrong results/error if files have the same paths, but are not the same.")
    parser.add_argument("-m", "--match_score", type=float, required=False, default=1, help="Match score. Default = 1. See biopyhton Align.PairwiseAlignment documentation.")
    parser.add_argument("-x", "--mismatch_score", type=float, required=False, default=0, help="Mismatach score. Default = 0. See biopyhton Align.PairwiseAlignment documentation.")
    parser.add_argument("-g", "--open_gap_score", type=float, required=False, default=0, help="Open gap score. Default = 0. See biopyhton Align.PairwiseAlignment documentation.")
    parser.add_argument("-e", "--extend_gap_score", type=float, required=False, default=0, help="Extended gap score. Default = 0. See biopyhton Align.PairwiseAlignment documentation.")
    parser.add_argument("-v", "--verbose", action="store_true", required=False, default=False, help="If given, more information is printed to the console")
    parser.add_argument("-c", "--cores", type=int, required=False, default=1, help="Number of cores used for the analysis. Default = 1")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = get_arguments()
    Pipeline(args)
