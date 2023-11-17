import sqlite3
from os import path, remove
from typing import List
import time

import numpy as np
import tqdm

def worker(line):
    if not line[0] == "#":
        # Seperate fields. Attribute endings get truncated, in order to use " as seperator
        chrom, _, feature, start, end, _, strand, _, attr = line.strip().split("\t")

        #  Check if the end of the attribute is either ";" or ". If not, raise an error
        if attr[-1] == ";":
            attr = attr[:-1]
        elif attr[-1] == '"':
            pass
        else:
            raise ValueError(
                'Was expecting the attribute ending to be either ";" or ". Check Ending of the attributes'
            )
        
        # Split the attributes into key-value pairs
        attr = attr.replace('"', "").replace("; ", ";").split(";")
        attr_dict = {}
        for pair in attr:
            key, value = pair.split(" ", maxsplit=1)
            attr_dict[key] = value
        attr = attr_dict

        # If the feature is a transcript and the transcript should be included in the database insert transcript data into the database
        if "transcript_id" in attr.keys():
            if included_transcripts is None or attr["transcript_id"] in included_transcripts:
                if feature == "transcript":
                
                    connection.execute(
                        f"INSERT INTO transcripts (transcript_id, chromosome, strand, start, end, sequence) VALUES (?, ?, ?, ?, ?, ?);",
                        (
                            attr["transcript_id"], 
                            chrom, 
                            strand, 
                            start, 
                            end, 
                            "",
                        ),
                    )

                # If the feature is an exon and the transcript should be included in the database insert exon data into the database
                elif feature == "exon":
                    connection.execute(
                        f"INSERT INTO exons (transcript_id, chromosome, strand, start, end) VALUES (?, ?, ?, ?, ?);",
                        (
                            attr["transcript_id"],
                            chrom,
                            strand,
                            start,
                            end,
                        ),
                    )

                # Commit changes to the database and close the connection
                connection.commit()

class AnnotationDB:
    """
    Description:
        A class to handle the creation and management of a SQLite database for GTF annotation files.

    Methods:
        create_annotation_db(gtf_path: str, tempdir: str, included_transcripts: np.ndarray) -> str:
            Creates a SQLite database for the given GTF file. The database will be created in the specified tempdir.
            The database will contain two tables: transcripts and exons.
        parse_gtf(gtf_path: str, db_path: str, included_transcripts: np.ndarray) -> None:
            Parses a gtf file and inserts the data into the database.
    """

    def create_annotation_db(gtf_path: str, tempdir: str, included_transcripts: np.ndarray, args) -> str:
        """
        Description:
            Creates a SQLite database for the given GTF file. The database will be created in the specified tempdir.
            The database will contain two tables: transcripts and exons.

            The transcripts table contains the following columns:
                transcript_id: TEXT
                chromosome: TEXT
                strand: TEXT
                start: INTEGER
                end: INTEGER
                sequence: TEXT
                
            The exons table contains the following columns:
                transcript_id: TEXT
                chromosome: TEXT
                strand: TEXT
                start: INTEGER
                end: INTEGER
                exon_number: INTEGER


        Parameters:
            gtf_path (str):
                Path to the GTF file.
            tempdir (str):
                Path to the tempdir where the database should be created in.
            included_transcripts (np.ndarray):
                A boolean array indicating which transcripts should be included in the database.
            args (argparse.ArgumentParser object):
                The parsed command-line arguments

        Returns:
            str
                Path to the created database.
        """
        # Get the name of the annotation from the gtf file path
        annotaion_name = gtf_path.split("/")[-1].split(".")[0]

        # Create the path to the database file
        db_path = tempdir + "/" + annotaion_name + ".db"

        # If the database already exists, return its path
        if path.exists(db_path):
            return db_path
        
        # Otherwise, create the database and its tables
        else:

            # Create the tables
            with sqlite3.connect(db_path) as connection:
                cursor = connection.cursor()
                cursor.execute(
                    "CREATE TABLE IF NOT EXISTS transcripts (transcript_id text PRIMARY KEY, chromosome text, strand text, start integer, end integer, sequence text);")
                cursor.execute(
                    "CREATE TABLE IF NOT EXISTS exons (transcript_id text, chromosome text, strand text, start integer, end integer);")
                
            # Parse the GTF file and insert the data into the database
            AnnotationDB.parse_gtf(gtf_path, db_path, included_transcripts, args)
            connection.commit()
            connection.close()

            return db_path

    def parse_gtf(gtf_path: str, db_path: str, included_transcripts: np.ndarray, args) -> None:
        """
        Description:
            Parses a gtf file and inserts the data into the database.
            The attribute endings are expected to be either ";" or ".

        Parameters:
            gtf_path (str):
                Path to the gtf file
            db_path (str):
                Path to the database
            included_transcripts (np.ndarray):
                A boolean array indicating which transcripts should be included in the database.
            args (argparse.ArgumentParser object):
                The parsed command-line arguments

        Returns:
            None

        Raises:
            ValueError
                If the attribute ending is not ";" or ".
        """
        # Open a connection to the database
        with sqlite3.connect(db_path) as connection:

            with open(gtf_path, "r") as gtf:
                # Iterate over lines in gtf
                # If verbose is True, use tqdm to print progress information
                lines = gtf.readlines()
                if args.verbose:
                    lines = tqdm.tqdm(lines, leave=False, desc="GTF lines")
                else:
                    lines = lines

                for line in lines:
                    if not line[0] == "#":
                        # Seperate fields. Attribute endings get truncated, in order to use " as seperator
                        chrom, _, feature, start, end, _, strand, _, attr = line.strip().split("\t")

                        #  Check if the end of the attribute is either ";" or ". If not, raise an error
                        if attr[-1] == ";":
                            attr = attr[:-1]
                        elif attr[-1] == '"':
                            pass
                        else:
                            raise ValueError(
                                'Was expecting the attribute ending to be either ";" or ". Check Ending of the attributes'
                            )
                        
                        # Split the attributes into key-value pairs
                        attr = attr.replace('"', "").replace("; ", ";").split(";")
                        attr_dict = {}
                        for pair in attr:
                            key, value = pair.split(" ", maxsplit=1)
                            attr_dict[key] = value
                        attr = attr_dict

                        # If the feature is a transcript and the transcript should be included in the database insert transcript data into the database
                        if "transcript_id" in attr.keys():
                            if included_transcripts is None or attr["transcript_id"] in included_transcripts:
                                if feature == "transcript":
                                
                                    connection.execute(
                                        f"INSERT INTO transcripts (transcript_id, chromosome, strand, start, end, sequence) VALUES (?, ?, ?, ?, ?, ?);",
                                        (
                                            attr["transcript_id"], 
                                            chrom, 
                                            strand, 
                                            start, 
                                            end, 
                                            "",
                                        ),
                                    )

                                # If the feature is an exon and the transcript should be included in the database insert exon data into the database
                                elif feature == "exon":
                                    connection.execute(
                                        f"INSERT INTO exons (transcript_id, chromosome, strand, start, end) VALUES (?, ?, ?, ?, ?);",
                                        (
                                            attr["transcript_id"],
                                            chrom,
                                            strand,
                                            start,
                                            end,
                                        ),
                                    )

                                # Commit changes to the database and close the connection
                                connection.commit()



class AnalysisDB:
    """
    Description:
        A class to handle the creation and management of a SQLite database for analysis of annotation files.

    Methods:
        create_analysis_db(outdir: str, reference_annotation_db_path: str, query_annotation_db_paths: List[str]) -> str:
            Creates an analysis database in the specified directory
    """

    def create_analysis_db(outdir: str, reference_annotation_db_path: str, query_annotation_db_paths: List[str]) -> str:
        """
        Creates an analysis database in the specified directory with the tables:
            max_similarity:
                transcript_id: TEXT
                query1: INTEGER
                query2: INTEGER
                ...
            query1_transcripts:
                transcript_id: TEXT
                max_similarity: REAL
            query1_results:
                threshold: REAL
                tp: INTEGER
                fp: INTEGER
                fn: INTEGER
                fdr: REAL
                tpr: REAL
            query2_transcripts:
                transcript_id: TEXT
                max_similarity: REAL
            query2_results:
                threshold: REAL
                tp: INTEGER
                fp: INTEGER
                fn: INTEGER
                fdr: REAL
                tpr: REAL
            ...

        Parameters:
            outdir (str):
                The path to the directory where the analysis database should be created.
            reference_annotation_db_path (str):
                The path to the reference annotation database.
            query_annotation_db_paths (List[str]):
                A list of paths to the query annotation files.

        Returns:
            str
                The path to the analysis database.
        """
        # Connect to the reference annotation database
        with sqlite3.connect(reference_annotation_db_path) as reference_annotation_db:

            # Connect to the query annotation databases
            query_annotation_dbs = [sqlite3.connect(query_annotation_db_path) for query_annotation_db_path in query_annotation_db_paths]

            # Get the names of the query annotation databases
            query_annotation_names = [query_annotation_db_path.split("/")[-1].split(".")[0] for query_annotation_db_path in query_annotation_db_paths]

            # Create the path to the analysis database
            db_path = path.join(outdir, "analysis.db")

            # Remove the analysis database if it already exists
            if path.exists(db_path):
                remove(db_path)

            # Connect to the analysis database
            with sqlite3.connect(db_path) as connection:

                # Create the max_similarity table
                max_similarity_columns = ["transcript_id"] + query_annotation_names
                max_similarity_query = f"CREATE TABLE IF NOT EXISTS max_similarity ({', '.join(max_similarity_columns)});"
                connection.execute(max_similarity_query)

                # Get the transcript ids from the reference annotation database
                transcript_ids = (reference_annotation_db.cursor().execute("SELECT transcript_id FROM transcripts").fetchall())

                # Insert the transcript ids into the max_similarity table
                values = [[transcript_id[0]] + [0] * len(query_annotation_names) for transcript_id in transcript_ids]

                connection.executemany(
                    f"INSERT INTO max_similarity ({', '.join(max_similarity_columns)}) VALUES ({', '.join(['?'] * len(max_similarity_columns))})",
                    values,
                )

                # Create the query transcripts and results tables
                for query_connection, query_name in zip(query_annotation_dbs, query_annotation_names):
                    
                    # Create the query transcripts table
                    connection.execute(
                        f"CREATE TABLE {query_name}_transcripts (transcript_id TEXT UNIQUE, max_similarity REAL)")
                    
                    # Insert the transcript ids into the query transcripts table
                    transcripts = query_connection.execute(
                        f"SELECT transcript_id FROM transcripts").fetchall()
                    
                    # Insert the transcript ids into the query transcripts table
                    values = [(transcript[0], 0) for transcript in transcripts]
                    
                    # Insert the transcript ids into the query transcripts table
                    connection.cursor().executemany(
                        f"INSERT OR IGNORE INTO {query_name}_transcripts (transcript_id, max_similarity) VALUES (?, ?)",
                        values,
                    )

                    # Create the query results table
                    connection.execute(
                        f"CREATE TABLE IF NOT EXISTS {query_name}_results (threshold REAL, tp_fdr INTEGER, tp_tpr INTEGER, fp INTEGER, fn INTEGER, fdr REAL, tpr REAL);")

                    # Commit changes to the query annotation database
                    query_connection.commit()
                    query_connection.close()

                # Commit changes to the analysis database
                connection.commit()

        # Return the path to the analysis database
        return db_path
