import sqlite3
from concurrent.futures import ProcessPoolExecutor
from time import process_time
from typing import Dict, List

import pandas as pd
import tqdm

from src.sequence_handler import Sequences


def worker(aruments):
    """
    Description:
        This function takes in a tuple of arguments and performs the following tasks:
        1. Connects to three SQLite databases: reference_connection, query_connection, and analysis_db_connection.
        2. Finds possible matches between a query and a subset of reference transcripts.
        3. Compares the query sequence with each possible match and calculates the similarity score.
        4. Updates the analysis_db_connection with the maximum similarity score for each transcript.
        5. Updates the analysis_db_connection with the maximum similarity score for the query.
        6. Closes the database connections.

    Parameters:
        aruments (tuple): 
            A tuple of arguments containing the following:

                reference_db_path (str): 
                    The path to the reference SQLite database.
                query_db_path (str): 
                    The path to the query SQLite database.
                analysis_db_path (str): 
                    The path to the analysis SQLite database.
                subset_reference (list): 
                    A list of reference transcripts to search for possible matches.
                chromosome_sequences (dict): 
                    A dictionary containing chromosome sequences.
                query_name (str): 
                    The name of the query.
                alignment_args (tuple): 
                    A tuple containing alignment arguments.
                query (tuple): 
                    A tuple containing information about the query transcript.

    Returns:
        None
    """

    # unpacking the arguments
    (
        reference_db_path,
        query_db_path,
        analysis_db_path,
        subset_reference,
        chromosome_sequences,
        query_name,
        alignment_args,
        query,
    ) = aruments

    # connecting to the databases
    with (
        sqlite3.connect(reference_db_path, timeout=1) as reference_connection,
        sqlite3.connect(analysis_db_path, timeout=1) as analysis_db_connection
    ):


        # finding possible matches between a query and a subset of reference transcripts
        if query[2] == "+":
            possible_matches = [reference for reference in subset_reference if (
                    ((reference[3] <= query[3]) and (query[3] <= reference[4]))
                    or ((query[3] <= reference[3]) and (reference[3] <= query[4]))
                )
            ]
        else:
            possible_matches = [reference for reference in subset_reference if (
                    ((reference[3] >= query[3]) and (query[3] >= reference[4]))
                    or ((query[3] >= reference[3]) and (reference[3] >= query[4]))
                )
            ]

        # getting the query sequence
        query_sequence = Sequences.get_sequence(query[0], query_db_path, chromosome_sequences)

        # updating the analysis_db_connection with the maximum similarity score for each transcript
        # if there are no possible matches, the similarity score is set to 0
        if len(possible_matches) == 0:
            similarity = 0
            start_time = process_time()
            while True:
                try:
                    analysis_db_connection.execute(
                        f"INSERT INTO {query_name}_transcripts (transcript_id, max_similarity) VALUES (?, ?)", (query[0], similarity))
                    analysis_db_connection.commit()
                except sqlite3.Error:
                    if process_time() - start_time > 3:
                        raise sqlite3.OperationalError("Timeout due to database lock")
                    else:
                        break
        
        # if there are possible matches, the similarity score is calculated and updated in the analysis_db_connection
        else:

            # loop through each possible match
            for match in possible_matches:
                start_time = process_time()

                # try to get the reference sequence from the database
                try:
                    reference_sequence = reference_connection.execute(
                        "SELECT sequence FROM transcripts WHERE transcript_id = ? ",(match[0],)).fetchone()[0]
                except sqlite3.Error:
                    break

                # if the reference sequence is not in the database, get it from the chromosome sequence and update the database
                if not reference_sequence:
                    reference_sequence = Sequences.get_sequence(match[0], reference_db_path, chromosome_sequences)
                    
                    # update the database with the reference sequence
                    start_time = process_time()
                    while True:
                        try:
                            reference_connection.execute(
                                "UPDATE transcripts SET sequence = ? WHERE transcript_id = ? ", (str(reference_sequence), match[0]))
                            reference_connection.commit()
                        except sqlite3.Error as e:
                            if process_time() - start_time > 3:
                                raise sqlite3.OperationalError(f"Timeout due to database lock {e}")
                        else:
                            break

            # calculate the similarity score
            similarity = Sequences.compare_sequences(reference_sequence, query_sequence, alignment_args)
            
            # update the analysis_db_connection with the maximum similarity score for each transcript
            start_time = process_time()
            while True:
                try:
                    analysis_db_connection.execute(
                        f"""
                        UPDATE max_similarity
                        SET {query_name} = ?
                        WHERE transcript_id = ?
                        AND ? > (
                            SELECT {query_name}
                            FROM max_similarity
                            WHERE transcript_id = ?
                        )
                    """,(similarity, match[0], similarity, match[0]))
                    analysis_db_connection.commit()
                except sqlite3.Error as e:
                    if process_time() - start_time > 3:
                        raise sqlite3.OperationalError(
                            "Timeout due to database lock" + str(e)
                        )
                else:
                    break

            # update the analysis_db_connection with the maximum similarity score for the query
            start_time = process_time()
            while True:
                try:
                    analysis_db_connection.execute(
                        f"""
                        UPDATE {query_name}_transcripts
                        SET max_similarity = ?
                        WHERE transcript_id = ?
                        AND max_similarity < ?
                        """,(similarity, query[0], similarity))
                    analysis_db_connection.commit()

                except sqlite3.Error as e:
                    if process_time() - start_time > 3:
                        raise sqlite3.OperationalError(
                            "Timeout due to database lock" + str(e)
                        )
                else:
                    break



class Analysis:
    """
    Description:
        A class for comparing reference annotation with a query annotation and calculating metrics for the analysis.

    Methods:
        compare_annotations(reference_db_path: str, query_db_path: str, analysis_db_path: str, chromosome_sequences: Dict[str, str], verbose: bool, args,) -> None:
            Compares the reference annotation with a query annotation.
        calculate_metrics(analysis_db_path: str, query_annotations: List[str]) -> None:
            Calculates the metrics for the analysis and fills the query_name_results tables For Threshold from 0.75 to 1.0 in steps of 0.01 TPR and FDR are calculated.
        produce_csv(analysis_db_path: str, query_annotation_names: List[str], output_dir: str) -> pd.DataFrame:
            Generate a CSV file containing analysis results for the given query annotation names and a combined CSV file containing the results for all query annotations.
    """

    def compare_annotations(reference_db_path: str, query_db_path: str, analysis_db_path: str, chromosome_sequences: Dict[str, str], args,) -> None:
        """
        Compares the reference annotation with a query annotation.

        Parameters:
            reference_db_path (str):
                The path to the reference annotation database.
            query_db_path (str):
                The path to the query annotation database.
            analysis_db_path (str):
                The path to the analysis database.
            chromosome_sequences (dict):
                A dictionary mapping chromosome names to their sequences.
            args (argparse.ArgumentParser object):
                An argparse.ArgumentParser object containing the command line arguments.

        Returns:
            None
        """
        verbose = args.verbose

        # Connect to the reference and query databases
        with (
            sqlite3.connect(reference_db_path, timeout=1) as reference_connection_high_level,
            sqlite3.connect(query_db_path, timeout=1) as query_connection_high_level,
        ):

            # Get the name of the query database
            query_name = query_db_path.split("/")[-1].split(".")[0]

            # Create a tuple of alignment arguments
            alignment_args = (
                args.match_score,
                args.mismatch_score,
                args.open_gap_score,
                args.extend_gap_score,
            )

            # Iterate over chromosomes
            # If verbose is True, use tqdm to print progress information
            if verbose:
                chromosomes = tqdm.tqdm(
                    [
                        row[0]
                        for row in query_connection_high_level.execute(
                            "SELECT DISTINCT chromosome FROM transcripts").fetchall()
                    ],
                    leave=False,
                    desc="Chromosome",
                )
            else:
                chromosomes = [
                    row[0]
                    for row in query_connection_high_level.execute(
                        "SELECT DISTINCT chromosome FROM transcripts").fetchall()
                ]

            # Iterate over chromosomes
            for chromosome in chromosomes:
                chromosomes.set_postfix({"current chromosome": chromosome}) if verbose else None

                # Iterate over strands
                if verbose:
                    strands = tqdm.tqdm(["+", "-"], leave=False, desc="Strand")
                else:
                    strands = ["+", "-"]
                for strand in strands:
                    strands.set_postfix({"current_strand": strand}) if verbose else None

                    # Get the subset of transcripts for the current chromosome and strand from the query and reference databases
                    subset_query = query_connection_high_level.execute(
                        f"SELECT * FROM transcripts WHERE chromosome = ? AND strand = ?", (chromosome, strand)).fetchall()
                    subset_reference = reference_connection_high_level.execute(
                        f"SELECT * FROM transcripts WHERE chromosome = ? AND strand = ?",(chromosome, strand),).fetchall()
                    reference_connection_high_level.commit()
                    query_connection_high_level.commit()

                    # Use a process pool to compare each query transcript with the reference transcripts
                    with ProcessPoolExecutor(max_workers=args.cores) as excecuter:

                        # Create a list of tasks for the process pool
                        tasks = [
                            (
                                reference_db_path,
                                query_db_path,
                                analysis_db_path,
                                subset_reference,
                                chromosome_sequences,
                                query_name,
                                alignment_args,
                                query,
                            )
                            for query in subset_query
                        ]

                        # Use a process pool to compare each query transcript with the reference transcripts
                        # If verbose is True, use tqdm to print progress information
                        if verbose:
                            list(
                                tqdm.tqdm(
                                    excecuter.map(worker, tasks),
                                    total=len(tasks),
                                    leave=False,
                                    desc="Reads",
                                )
                            )

                        # If verbose is False, just use a process pool to compare each query transcript with the reference transcripts
                        else:
                            list(excecuter.map(worker, tasks))


    def calculate_metrics(analysis_db_path: str, query_annotations: List[str]) -> None:
        """
        Description:
            Calculates the metrics for the analysis and fills the query_name_results tables For Threshold from 0.75 to 1.0 in steps of 0.01 TPR and FDR are calculated.

        Definitions:
            TP_FDR: Number of true positives
                number of max_similarity values >= threshold in the max_similarity column of the respective query_name_transcripts table
            FP: Number of false positives
                number of max_similarity values < threshold in the max_similarity column of the respective query_name_transcripts table
            FDR: False discovery rate (How many of the predicted transcripts are not in the reference annotation)
                FP / (FP + TP)

            TP_TPR: Number of true positives
                number of max_similarity values >= threshold in respective query column of the max_similarity table
            FN: Number of false negatives
                number of max_similarity values < threshold in the respective query column of the max_similarity table
            TPR: True positive rate / Sensitivity (How many of the transcripts in the reference annotation are predicted)
                TP / (TP + FN)

        Parameters:
            analysis_db_path (str):
                The path to the analysis database.
            query_annotations (List[str]):
                A list of query annotations.

        Returns:
            None
        """

        # split the query_annotations to get the query names
        query_names = [query_annotation.split("/")[-1].split(".")[0] for query_annotation in query_annotations]

        # connect to the analysis database
        analysis_db_connection = sqlite3.connect(analysis_db_path)

        # loop through each query name
        for query_name in query_names:

            # loop through each threshold value from 1.0 to 0.75 in steps of 0.01
            for threshold in list(range(75, 101, 1))[::-1]:
                threshold_value = threshold / 100

                # calculate TP_FDR
                tp_fdr = analysis_db_connection.execute(
                    f"SELECT COUNT(*) FROM {query_name}_transcripts WHERE max_similarity >= ?",(threshold_value,),).fetchone()[0]
                
                # calculate FP
                fp = analysis_db_connection.execute(
                    f"SELECT COUNT(*) FROM {query_name}_transcripts WHERE max_similarity < ?",(threshold_value,),).fetchone()[0]
                # calculate FDR
                fdr = fp / (fp + tp_fdr)

                # calculate TP_TPR
                tp_tpr = analysis_db_connection.execute(
                    f"SELECT COUNT(*) FROM max_similarity WHERE {query_name} >= ?",
                    (threshold_value,),
                ).fetchone()[0]

                # calculate FN
                fn = analysis_db_connection.execute(
                    f"SELECT COUNT(*) FROM max_similarity WHERE {query_name} < ?",
                    (threshold_value,),
                ).fetchone()[0]

                # calculate TPR
                tpr = tp_tpr / (tp_tpr + fn)

                # insert the calculated values into the query_name_results table
                analysis_db_connection.execute(
                    f"INSERT INTO {query_name}_results (threshold, tp_fdr, tp_tpr, fp, fn, tpr, fdr) VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (threshold_value, tp_fdr, tp_tpr, fp, fn, tpr, fdr),
                )

        # commit the changes to the database and close the connection
        analysis_db_connection.commit()
        analysis_db_connection.close()

    def produce_csv(analysis_db_path: str, query_annotation_names: List[str], output_dir: str) -> pd.DataFrame:
        """
        Description:
            Generate a CSV file containing analysis results for the given query annotation names and a combined
            CSV file containing the results for all query annotations.

        Parameters: 
            analysis_db_path (str):
                The path to the analysis database.
            query_annotation_names (List[str]):
                A list of query annotation names to include in the CSV.
            output_dir (str):
                The directory to write the CSV files to.

        Returns:
            pd.DataFrame
                A pandas DataFrame containing the combined results for all query annotation names.
        """

        with sqlite3.connect(analysis_db_path) as analysis_db_connection:

            # Initialize variables
            query_annotation_results = {}
            thresholds = []

            # Loop through each query annotation name
            for query_annotation_name in query_annotation_names:

                # Retrieve query annotation results from the analysis database
                query_annotation_result = analysis_db_connection.execute(
                    f"SELECT * FROM {query_annotation_name}_results").fetchall()

                # Convert query annotation results to a pandas DataFrame
                df = pd.DataFrame(
                    query_annotation_result,
                    columns=["threshold", "tp_fdr", "tp_tpr", "fp", "fn", "tpr", "fdr"],
                )

                # Write query annotation results to a CSV file
                df.to_csv(f"{output_dir}/{query_annotation_name}_results.csv", index=False)

                # Store query annotation results in a dictionary
                query_annotation_results[query_annotation_name] = df[["threshold", "fdr", "tpr"]]

                # Store thresholds if they have not been stored yet
                if len(thresholds) == 0:
                    thresholds = df["threshold"].tolist()


        # Combine query annotation results into a single pandas DataFrame
        combined_df = pd.DataFrame({"threshold": thresholds})
        for query_annotation_name in query_annotation_names:
            if query_annotation_name in query_annotation_results:
                combined_df = pd.merge(
                    combined_df,
                    query_annotation_results[query_annotation_name],
                    on="threshold",
                    how="left",
                )
                combined_df = combined_df.rename(
                    columns={
                        "fdr": f"{query_annotation_name}_fdr",
                        "tpr": f"{query_annotation_name}_tpr",
                    }
                )
            else:
                combined_df[f"{query_annotation_name}_fdr"] = ""
                combined_df[f"{query_annotation_name}_tpr"] = ""

        # Write combined query annotation results to a CSV file
        combined_df.to_csv(f"{output_dir}/combined_results.csv", index=False)

        # Return the combined query annotation results as a pandas DataFrame
        return combined_df
