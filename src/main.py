import os
import shutil
from pathlib import Path
from shutil import rmtree
from time import process_time

import numpy as np
import tqdm

from src.analysis import Analysis
from src.database_handler import AnalysisDB, AnnotationDB
from src.sequence_handler import FastaParser


class Pipeline:
    """
    Description:
        Pipeline class that represents the main pipeline for the NoTraDAM tool.

    Attributes:
        temp_dir (str):
            Temporary directory for databases.
        figures_dir (str):
            Directory for figures.
        transcript_ids (numpy.ndarray or None):
            Array of transcript IDs that should be included. If None, all transcripts are included.
        query_annotation_db_paths (dict):
            Dictionary with annotation paths as keys and database connections as values.
        reference_annotation_db_path (str):
            Path to the reference annotation database.
        chromosome_sequences (dict):
            Dictionary with chromosome names as keys and sequences as values.
        analysis_db_path (str):

    Methods:
        __init__(self, args)
            Initializes the Pipeline object.
        create_output_dir(self, args)
            Creates output, temp and figures directories, if they do not already exist.
        parse_transcript_ids(self, args)
            Parses the transcript text file. If none is given, then transcript_ids is None.
        create_annotation_dbs(self, args)
            Creates annotation databases.
        get_chromosome_sequences(self, args)
            Calls FastaParser to parse the chromosome fasta file.
        create_analysis_db(self, args)
            Creates analysis database.
        compare_reference_with_queries(self, args)
            Compares the reference annotation with each of the query annotations provided in the arguments.
        calculate_statistics(self, args)
            Calculates statistics for the analysis database.
    """

    def __init__(self, args):
        """
        Description:
            Initializes the Pipeline object and calls the different methods in the pipeline.

        Parameters:
            args (argparse.ArgumentParser object):
                The parsed command-line arguments

        Returns:
            None
        """
        # Create output directory
        self.create_output_dir(args)

        # Parse transcript ids
        print(
            "\n\n(1/5) Parsing transcript ids -------------------------------------\n"
        ) if args.verbose else None
        self.parse_transcript_ids(args)

        # Parse annotations
        print(
            "\n(2/5) Parsing annotations ----------------------------------------\n"
        ) if args.verbose else None
        self.create_annotation_dbs(args)
        self.create_analysis_db(args)

        # Parse chromosome sequences
        print(
            "\n(3/5) Parsing chromosome sequences -------------------------------\n"
        ) if args.verbose else None
        self.get_chromosome_sequences(args)

        # Compare reference with queries
        print(
            "\n(4/5) Comparing reference with queries ---------------------------\n"
        ) if args.verbose else None
        self.compare_reference_with_queries(args)

        # Calculate statistics
        print(
            "\n(5/5) Calculating statistics -------------------------------------\n"
        ) if args.verbose else None
        self.calculate_statistics(args)

        # Remove temporary directory if keep_temporary is False
        if not args.keep_temporary:
            rmtree(self.temp_dir)

    def create_output_dir(self, args):
        """
        Description:
            Create output, temp and figures directories, if it not already exists.

        Parameters:
        args (argparse.ArgumentParser object):
            The parsed command-line arguments

        Returns:
            None
        """
        # If output_directory is not provided, set it to current working directory + "/NoTraDAM_output"
        if args.output_directory is None:
            args.output_directory = os.getcwd()
            output_directory = f"{args.output_directory}/NoTraDAM_output"
        else:
            output_directory = args.output_directory


        self.temp_dir = f"{output_directory}/temp"

        # Create output_directory and self.temp_dir if they don't already exist
        Path(output_directory).mkdir(parents=True, exist_ok=True)
        Path(self.temp_dir).mkdir(exist_ok=True)

    def parse_transcript_ids(self, args):
        """
        Description:
            Parses transcirpt text file. If none is given, then transcript_ids is none.

        Parameters:
            args (argparse.ArgumentParser object):
                The parsed command-line arguments

        Returns:
            None
        """
        self.transcript_ids = {}

        # If no transcript_ids file is provided, set the reference annotation and query annotations to None
        if args.transcript_ids is None:
            self.transcript_ids[self.reference_annotation] = None
            for annotation in args.query_annotations:
                self.transcript_ids[annotation] = None

        # If a transcript_ids file is provided, parse it
        else:
            with open(args.transcript_ids, "r") as f:
                current_path = None
                current_transcripts = []
                
                for line in f:
                    line = line.strip()
                    if line.startswith("@"):
                        if current_path is not None:
                            self.transcript_ids[current_path] = np.array(current_transcripts)
                        current_path = line[1:]
                        current_transcripts = []
                    elif line:
                        current_transcripts.append(line)
                        
                if current_path is not None:
                    self.transcript_ids[current_path] = np.array(current_transcripts)


    def create_annotation_dbs(self, args):
        """
        Description:
            Creates annotation databases.

        Parameters:
            args (argparse.ArgumentParser object):
                The parsed command-line arguments

        Returns:
            None
        """
        # Start the timer if verbose mode is on
        start = process_time() if args.verbose else None

        # Initialize the query_annotation_db_paths dictionary
        self.query_annotation_db_paths = {}

        # Create the annotation databases
        if args.verbose:
            annotations = tqdm.tqdm([args.reference_annotation]+args.query_annotations, leave=False, desc="Query Annotation")
        else:
            annotations = [args.reference_annotation]+args.query_annotations

        for i, annotation in enumerate(annotations):
            annotations.set_postfix({"current query annotation": annotation}) if args.verbose else None
            if i == 0:
                self.reference_annotation_db_path = AnnotationDB.create_annotation_db(
                annotation, self.temp_dir, self.transcript_ids[annotation], args)
            else:
                self.query_annotation_db_paths[annotation] = AnnotationDB.create_annotation_db(
                    annotation, self.temp_dir, self.transcript_ids[annotation], args)

        # Print the number of query annotations and processing time if verbose mode is on
        print(
            f"\tNumber of query annotations: {len(self.query_annotation_db_paths)}\n\tprocessing time: {(process_time()-start):0.1e} s"
        ) if args.verbose else None

    def get_chromosome_sequences(self, args):
        """
        Description:
            Calling FastaParser to parse the chromosome fasta file.

        Parameters:
            args (argparse.ArgumentParser object):
                The parsed command-line arguments

        Returns:
            None
        """
        # Start the timer if verbose mode is on
        start = process_time() if args.verbose else None

        # Parse the chromosome fasta file using FastaParser
        self.chromosome_sequences = FastaParser.parse_fasta(args.genome_fasta)

        # Print the number of chromosomes and processing time if verbose mode is on
        print(
            f"\tNumber of chromosomes: {len(self.chromosome_sequences)}\n\tprocessing time: {(process_time()-start):0.1e} s"
        ) if args.verbose else None

    def create_analysis_db(self, args):
        """
        Description:
            Creates analysis database.

        Parameters:
            args (argparse.ArgumentParser object):
                The parsed command-line arguments

        Returns:
            None
        """
        # Set the analysis_db_path attribute to the path of the newly created analysis database
        self.analysis_db_path = AnalysisDB.create_analysis_db(
            args.output_directory,
            self.reference_annotation_db_path,
            self.query_annotation_db_paths.values(),
        )

    def compare_reference_with_queries(self, args):
        """
        Description:
            Compares the reference annotation with each of the query annotations provided in the arguments.

        Parameters:
            args (argparse.ArgumentParser object):
                The parsed command-line arguments

        Returns:
            None
        """
        # Start the timer if verbose mode is on
        start = process_time() if args.verbose else None

        # Create a progress bar for the query annotations
        query_annotations = tqdm.tqdm(args.query_annotations, leave=False, desc="Annotations") if args.verbose else args.query_annotations
        # Loop through each query annotation and compare it with the reference annotation
        for query_annotation in query_annotations:
            # Update the progress bar with the current query strand
            if args.verbose:
                query_annotations.set_postfix({"current_strand": query_annotation}) 

            # Compare the annotations
            Analysis.compare_annotations(
                self.reference_annotation_db_path,
                self.query_annotation_db_paths[query_annotation],
                self.analysis_db_path,
                self.chromosome_sequences,
                args,
            )

        # Print the processing time if verbose mode is on
        print(
            f"\tprocessing time: {(process_time()-start):0.1e} s"
        ) if args.verbose else None

    def calculate_statistics(self, args):
        """
        Description:
            Calculates statistics for the analysis database.

        Parameters:
            args (argparse.ArgumentParser object):
                The parsed command-line arguments

        Returns:
            None
        """
        # Start the timer if verbose mode is on
        start = process_time() if args.verbose else None

        # Extract the annotation names from the query annotations
        query_annotation_names = [annotation.split("/")[-1].split(".")[0] for annotation in args.query_annotations]

        # Calculate the metrics for the analysis database
        Analysis.calculate_metrics(self.analysis_db_path, query_annotation_names)

        # Produce a CSV file with the analysis results
        df = Analysis.produce_csv(self.analysis_db_path, query_annotation_names, args.output_directory)

        # Print the results and processing time if verbose mode is on
        print(
            f"\tResults are stored in {args.output_directory}\n\tprocessing time: {(process_time()-start):0.1e} s\n\n\n\n{df.head().to_markdown(index=False)}"
        ) if args.verbose else None
