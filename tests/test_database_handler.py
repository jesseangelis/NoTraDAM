import sqlite3
from os import path, remove
from shutil import rmtree

import numpy as np
import pytest

from src.database_handler import AnalysisDB, AnnotationDB


class TestAnnotations:
    @pytest.fixture
    def generate_gtf_files(self):
        reference_annotation_content = """##description:
##provider: ...
##contact: ...
##format: gtf
##date: 2023-01-01
chr1	TEST	gene	1	3	.	+	.	gene_id "ENSG11"; gene_type "lncRNA"; gene_name "Gene11";
chr1	TEST	transcript	1	3	.	+	.	gene_id "ENSG11"; transcript_id "ENST111"; gene_type "lncRNA"; gene_name "Gene11"; transcript_type "lncRNA"; transcript_name "Transcript111";
chr1	TEST	exon	1	3	.	+	.	gene_id "ENSG11"; transcript_id "ENST111"; gene_type "lncRNA"; gene_name "Gene1"; transcript_type "lncRNA"; transcript_name "Transcript111"; exon_number 1; exon_id "Exon1111";
chr1	TEST	gene	200	205	.	+	.	gene_id "ENSG12"; gene_type "lncRNA"; gene_name "Gene12";
chr1	TEST	transcript	200	205	.	+	.	gene_id "ENSG12"; transcript_id "ENST121"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript121";
chr1	TEST	exon	200	201	.	+	.	gene_id "ENSG12"; transcript_id "ENST121"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript121"; exon_number 1; exon_id "Exon1211";
chr1	TEST	exon	204	205	.	+	.	gene_id "ENSG12"; transcript_id "ENST121"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript121"; exon_number 2; exon_id "Exon1213";
chr1	TEST	transcript	200	203	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122";
chr1	TEST	exon	200	201	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122"; exon_number 1; exon_id "Exon1221";
chr1	TEST	exon	204	203	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122"; exon_number 2; exon_id "Exon1222";
chr2	TEST	gene	105	100	.	-	.	gene_id "ENSG21"; gene_type "lncRNA"; gene_name "Gene21";
chr2	TEST	transcript	105	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211";
chr2	TEST	exon	105	104	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 1; exon_id "Exon2111";
chr2	TEST	exon	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 2; exon_id "Exon2112";
chr2	TEST	transcript	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript212";
chr2	TEST	exon	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 1; exon_id "Exon2112";
chr2	TEST	gene	300	302	.	+	.	gene_id "ENSG22"; gene_type "lncRNA"; gene_name "Gene22";
chr2	TEST	transcript	300	302	.	+	.	gene_id "ENSG22"; transcript_id "ENST221"; gene_type "lncRNA"; gene_name "Gene22"; transcript_type "lncRNA"; transcript_name "Transcript221";
chr2	TEST	exon	300	302	.	+	.	gene_id "ENSG22"; transcript_id "ENST221"; gene_type "lncRNA"; gene_name "Gene22"; transcript_type "lncRNA"; transcript_name "Transcript221"; exon_number 1; exon_id "Exon2211";
chr2	TEST	gene	302	300	.	-	.	gene_id "ENSG23"; gene_type "lncRNA"; gene_name "Gene23";
chr2	TEST	transcript	302	300	.	-	.	gene_id "ENSG23"; transcript_id "ENST231"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Transcript231";
chr2	TEST	exon	302	300	.	-	.	gene_id "ENSG23"; transcript_id "ENST231"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Transcript231"; exon_number 1; exon_id "Exon2311";"""
        with open("tests/reference_annotation.gtf", "w") as f:
            f.write(reference_annotation_content)
        yield
        if path.exists("tests/reference_annotation.gtf"):
            remove("tests/reference_annotation.gtf")

    @pytest.fixture
    def dummy_args_instance(self):
        class DummyArgs:
            def __init__(self):
                self.verbose = False
                self.keep_temporary = False
                self.match_score = 1
                self.mismatch_score = 0
                self.open_gap_score = 0
                self.extend_gap_score = 0
                self.cores = 1
                self.included_transcript_ids = None
                self.output_directory = "tests"
                self.reference_annotation = "tests/reference_annotation.gtf"
                self.query_annotations = [
                    "tests/query_annotation1.gtf",
                    "tests/query_annotation2.gtf",
                ]
                self.genome_fasta = "tests/genome.fasta"

        yield DummyArgs()
    
    def test_parse_gtf(self, generate_gtf_files, dummy_args_instance):
        """
        Description:
            Test case to check if the gtf is parsed correctly.

        Args:
            generate_gtf_files: A fixture that generates GTF files for testing.

        Returns:
            None
        """
        # Remove the test database if it exists
        if path.exists("tests/test.db"):
            remove("tests/test.db")
        # Create a connection to the test database
        connection = sqlite3.connect("tests/test.db")
        # Create the transcripts and exons tables
        connection.execute(
            "CREATE TABLE transcripts (transcript_id text, chromosome text, strand text, start integer, end integer, sequence text)"
        )
        connection.execute(
            "CREATE TABLE exons (transcript_id text, chromosome text, strand text, start integer, end integer)"
        )
        # Call the parse_gtf function
        AnnotationDB.parse_gtf(
            "tests/reference_annotation.gtf",
            "tests/test.db",
            np.array(
                ["ENST111", "ENST121", "ENST122", "ENST221", "ENST211", "ENST212"]
            ),
            dummy_args_instance
        )
        # Check if the data was inserted correctly into the transcripts table
        transcripts = connection.execute("SELECT * FROM transcripts").fetchall()
        assert transcripts == [
            ("ENST111", "chr1", "+", 1, 3, ""),
            ("ENST121", "chr1", "+", 200, 205, ""),
            ("ENST122", "chr1", "+", 200, 203, ""),
            ("ENST211", "chr2", "-", 105, 100, ""),
            ("ENST212", "chr2", "-", 102, 100, ""),
            ("ENST221", "chr2", "+", 300, 302, ""),
        ]
        # Check if the data was inserted correctly into the exons table
        exons = connection.execute("SELECT * FROM exons").fetchall()
        assert exons == [
            ("ENST111", "chr1", "+", 1, 3),
            ("ENST121", "chr1", "+", 200, 201),
            ("ENST121", "chr1", "+", 204, 205),
            ("ENST122", "chr1", "+", 200, 201),
            ("ENST122", "chr1", "+", 204, 203),
            ("ENST211", "chr2", "-", 105, 104),
            ("ENST211", "chr2", "-", 102, 100),
            ("ENST212", "chr2", "-", 102, 100),
            ("ENST221", "chr2", "+", 300, 302),
        ]
        # Close the connection to the test database
        connection.close()
        # Remove the test database and output directory if they exist
        if path.exists("tests/test.db"):
            remove("tests/test.db")
        if path.exists("tests/output"):
            rmtree("tests/output")

    def test_create_annotation_db(self, generate_gtf_files, dummy_args_instance):
        """
        Description:
            Test case to check if the annotation database is created successfully.

        Args:
            generate_gtf_files: A fixture that generates GTF files for testing.

        Returns:
            None
        """
        # Remove the old reference annotation database if it exists
        if path.exists("tests/reference_annotaiton.gtf"):
            remove("tests/reference_annotation.db")

        # Create the annotation database
        db_path = AnnotationDB.create_annotation_db(
            "tests/reference_annotation.gtf", "tests", None, dummy_args_instance
        )

        # Check if the database connection is successful
        connection = sqlite3.connect(db_path)
        assert connection is not None

        # Check if the database file exists
        assert path.exists("tests/reference_annotation.db")

        # Check if the tables are created successfully
        tables = connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table';"
        ).fetchall()
        assert tables == [("transcripts",), ("exons",)]

        # Close the database connection
        connection.close()

        # Remove the reference annotation database and output directory
        if path.exists("tests/reference_annotation.db"):
            remove("tests/reference_annotation.db")
        if path.exists("tests/output"):
            rmtree("tests/output")


class TestAnalysisDB:
    @pytest.fixture
    def generte_reference_gtf(self):
        if path.exists("tests/reference_annotation.gtf"):
            remove("tests/reference_annotation.gtf")
        reference_annotation_content = """##description: Not in GTF: ENST121; Not in transcript ids: ENST211 
##provider: ...
##contact: ...
##format: gtf
##date: 2023-01-01
chr1	TEST	gene	1	3	.	+	.	gene_id "ENSG11"; gene_type "lncRNA"; gene_name "Gene11";
chr1	TEST	transcript	1	3	.	+	.	gene_id "ENSG11"; transcript_id "ENST111"; gene_type "lncRNA"; gene_name "Gene11"; transcript_type "lncRNA"; transcript_name "Transcript111";
chr1	TEST	exon	1	3	.	+	.	gene_id "ENSG11"; transcript_id "ENST111"; gene_type "lncRNA"; gene_name "Gene1"; transcript_type "lncRNA"; transcript_name "Transcript111"; exon_number 1; exon_id "Exon1111";
chr1	TEST	gene	200	205	.	+	.	gene_id "ENSG12"; gene_type "lncRNA"; gene_name "Gene12";
chr1	TEST	transcript	200	203	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122";
chr1	TEST	exon	200	201	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122"; exon_number 1; exon_id "Exon1221";
chr1	TEST	exon	203	205	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122"; exon_number 2; exon_id "Exon1222";
chr2	TEST	gene	105	100	.	-	.	gene_id "ENSG21"; gene_type "lncRNA"; gene_name "Gene21";
chr2	TEST	transcript	105	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211";
chr2	TEST	exon	105	104	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 1; exon_id "Exon2111";
chr2	TEST	exon	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 2; exon_id "Exon2112";
chr2	TEST	transcript	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript212";
chr2	TEST	exon	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 1; exon_id "Exon2112";
chr2	TEST	gene	300	302	.	+	.	gene_id "ENSG22"; gene_type "lncRNA"; gene_name "Gene22";
chr2	TEST	transcript	300	302	.	+	.	gene_id "ENSG22"; transcript_id "ENST221"; gene_type "lncRNA"; gene_name "Gene22"; transcript_type "lncRNA"; transcript_name "Transcript221";
chr2	TEST	exon	300	302	.	+	.	gene_id "ENSG22"; transcript_id "ENST221"; gene_type "lncRNA"; gene_name "Gene22"; transcript_type "lncRNA"; transcript_name "Transcript221"; exon_number 1; exon_id "Exon2211";
chr2	TEST	gene	302	300	.	-	.	gene_id "ENSG23"; gene_type "lncRNA"; gene_name "Gene23";
chr2	TEST	transcript	302	300	.	-	.	gene_id "ENSG23"; transcript_id "ENST231"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Transcript231";
chr2	TEST	exon	302	300	.	-	.	gene_id "ENSG23"; transcript_id "ENST231"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Transcript231"; exon_number 1; exon_id "Exon2311";"""
        with open("tests/reference_annotation.gtf", "w") as f:
            f.write(reference_annotation_content)
        yield
        if path.exists("tests/reference_annotation.gtf"):
            remove("tests/reference_annotation.gtf")

    @pytest.fixture
    def generte_query1_gtf(self):
        if path.exists("tests/query_annotation1.gtf"):
            remove("tests/query_annotation1.gtf")
        query_annotation1_content = """##description: Not in GTF: ENST111, ENST211; Not in transcript ids: ENST211; Added to GTF: ENSTWrong
##provider: ...
##contact: ...
##format: gtf
##date: 2023-01-01
chr1	TEST	gene	200	205	.	+	.	gene_id "ENSG12"; gene_type "lncRNA"; gene_name "Gene12";
chr1	TEST	transcript	200	205	.	+	.	gene_id "ENSG12"; transcript_id "ENST121"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript121";
chr1	TEST	exon	200	201	.	+	.	gene_id "ENSG12"; transcript_id "ENST121"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript121"; exon_number 1; exon_id "Exon1211";
chr1	TEST	exon	204	205	.	+	.	gene_id "ENSG12"; transcript_id "ENST121"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript121"; exon_number 2; exon_id "Exon1213";
chr1	TEST	transcript	200	203	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122";
chr1	TEST	exon	200	201	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122"; exon_number 1; exon_id "Exon1221";
chr1	TEST	exon	203	205	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122"; exon_number 2; exon_id "Exon1222";
chr2	TEST	gene	105	100	.	-	.	gene_id "ENSG21"; gene_type "lncRNA"; gene_name "Gene21";
chr2	TEST	transcript	105	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211";
chr2	TEST	exon	105	104	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 1; exon_id "Exon2111";
chr2	TEST	exon	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 2; exon_id "Exon2112";
chr2	TEST	transcript	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript212";
chr2	TEST	exon	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 1; exon_id "Exon2112";
chr2	TEST	gene	302	300	.	-	.	gene_id "ENSG23"; gene_type "lncRNA"; gene_name "Gene23";
chr2	TEST	transcript	302	300	.	-	.	gene_id "ENSG23"; transcript_id "ENST231"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Transcript231";
chr2	TEST	exon	302	300	.	-	.	gene_id "ENSG23"; transcript_id "ENST231"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Transcript231"; exon_number 1; exon_id "Exon2311";
chr2	TEST	gene	200	204	.	-	.	gene_id "ENSGWrong"; gene_type "lncRNA"; gene_name "Wrong";
chr2	TEST	transcript	200	204	.	-	.	gene_id "ENSGWrong"; transcript_id "ENSTWrong"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Wrong";
chr2	TEST	exon	200	204	.	-	.	gene_id "ENSGWrong"; transcript_id "ENSTWrong"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Wrong"; exon_number 1; exon_id "EXONWrong";"""
        with open("tests/query_annotation1.gtf", "w") as f:
            f.write(query_annotation1_content)
        yield
        if path.exists("tests/query_annotation1.gtf"):
            remove("tests/query_annotation1.gtf")

    @pytest.fixture
    def generte_query2_gtf(self):
        if path.exists("tests/query_annotation2.gtf"):
            remove("tests/query_annotation2.gtf")
        query_annotation2_content = """##description: Not in GTF: ENST121; 
##provider: ...
##contact: ...
##format: gtf
##date: 2023-01-01
chr1	TEST	gene	1	3	.	+	.	gene_id "ENSG11"; gene_type "lncRNA"; gene_name "Gene11";
chr1	TEST	transcript	1	3	.	+	.	gene_id "ENSG11"; transcript_id "ENST111"; gene_type "lncRNA"; gene_name "Gene11"; transcript_type "lncRNA"; transcript_name "Transcript111";
chr1	TEST	exon	1	3	.	+	.	gene_id "ENSG11"; transcript_id "ENST111"; gene_type "lncRNA"; gene_name "Gene1"; transcript_type "lncRNA"; transcript_name "Transcript111"; exon_number 1; exon_id "Exon1111";
chr1	TEST	gene	200	205	.	+	.	gene_id "ENSG12"; gene_type "lncRNA"; gene_name "Gene12";
chr1	TEST	transcript	200	203	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122";
chr1	TEST	exon	200	201	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122"; exon_number 1; exon_id "Exon1221";
chr1	TEST	exon	203	205	.	+	.	gene_id "ENSG12"; transcript_id "ENST122"; gene_type "lncRNA"; gene_name "Gene12"; transcript_type "lncRNA"; transcript_name "Transcript122"; exon_number 2; exon_id "Exon1222";
chr2	TEST	gene	105	100	.	-	.	gene_id "ENSG21"; gene_type "lncRNA"; gene_name "Gene21";
chr2	TEST	transcript	105	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211";
chr2	TEST	exon	105	104	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 1; exon_id "Exon2111";
chr2	TEST	exon	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 2; exon_id "Exon2112";
chr2	TEST	transcript	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript212";
chr2	TEST	exon	102	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 1; exon_id "Exon2112";
chr2	TEST	gene	300	302	.	+	.	gene_id "ENSG22"; gene_type "lncRNA"; gene_name "Gene22";
chr2	TEST	transcript	300	302	.	+	.	gene_id "ENSG22"; transcript_id "ENST221"; gene_type "lncRNA"; gene_name "Gene22"; transcript_type "lncRNA"; transcript_name "Transcript221";
chr2	TEST	exon	300	302	.	+	.	gene_id "ENSG22"; transcript_id "ENST221"; gene_type "lncRNA"; gene_name "Gene22"; transcript_type "lncRNA"; transcript_name "Transcript221"; exon_number 1; exon_id "Exon2211";
chr2	TEST	gene	302	300	.	-	.	gene_id "ENSG23"; gene_type "lncRNA"; gene_name "Gene23";
chr2	TEST	transcript	302	300	.	-	.	gene_id "ENSG23"; transcript_id "ENST231"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Transcript231";
chr2	TEST	exon	302	300	.	-	.	gene_id "ENSG23"; transcript_id "ENST231"; gene_type "lncRNA"; gene_name "Gene23"; transcript_type "lncRNA"; transcript_name "Transcript231"; exon_number 1; exon_id "Exon2311";"""
        with open("tests/query_annotation2.gtf", "w") as f:
            f.write(query_annotation2_content)
        yield
        if path.exists("tests/query_annotation2.gtf"):
            remove("tests/query_annotation2.gtf")
        
    @pytest.fixture
    def dummy_args_instance(self):
        class DummyArgs:
            def __init__(self):
                self.verbose = False
                self.keep_temporary = False
                self.match_score = 1
                self.mismatch_score = 0
                self.open_gap_score = 0
                self.extend_gap_score = 0
                self.cores = 1
                self.included_transcript_ids = None
                self.output_directory = "tests"
                self.reference_annotation = "tests/reference_annotation.gtf"
                self.query_annotations = [
                    "tests/query_annotation1.gtf",
                    "tests/query_annotation2.gtf",
                ]
                self.genome_fasta = "tests/genome.fasta"

        yield DummyArgs()

    def test_create_analysis_db(
        self, generte_reference_gtf, generte_query1_gtf, generte_query2_gtf, dummy_args_instance):
        """
        Description:
            Test case to check if the AnalysisDB.create_analysis_db method creates the expected database.

        Args:
            generte_reference_gtf (fixture): 
                A fixture that generates a reference GTF file.
            generte_query1_gtf (fixture):
                A fixture that generates a query1 GTF file.
            generte_query2_gtf (fixture):
                A fixture that generates a query2 GTF file.
        """
        # Remove existing analysis.db file if it exists
        if path.exists("tests/analysis.db"):
            remove("tests/analysis.db")

        # Create analysis database
        db_path = AnalysisDB.create_analysis_db(
            "tests",
            AnnotationDB.create_annotation_db(
                "tests/reference_annotation.gtf", "tests", None, dummy_args_instance
            ),
            [
                AnnotationDB.create_annotation_db(
                    "tests/query_annotation1.gtf", "tests", None, dummy_args_instance
                ),
                AnnotationDB.create_annotation_db(
                    "tests/query_annotation2.gtf", "tests", None, dummy_args_instance
                ),
            ],
        )

        # Check if analysis.db file was created
        assert path.exists("tests/analysis.db")

        # Check if expected tables were created in the database
        connection = sqlite3.connect(db_path)
        tables = connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table';"
        ).fetchall()
        expected_tables = [
            ("max_similarity",),
            ("query_annotation1_transcripts",),
            ("query_annotation1_results",),
            ("query_annotation2_transcripts",),
            ("query_annotation2_results",),
        ]
        assert tables == expected_tables, tables

        # Check if expected transcript ids were inserted into the max_similarity table
        transcript_ids = connection.execute(
            "SELECT transcript_id FROM max_similarity"
        ).fetchall()
        expected_ids = [
            ("ENST111"),
            ("ENST122"),
            ("ENST211"),
            ("ENST212"),
            ("ENST221"),
            ("ENST231"),
        ]
        for r, e in zip(transcript_ids, expected_ids):
            assert r[0] == e
        connection.close()

        # Remove created files
        if path.exists("tests/analysis.db"):
            remove("tests/analysis.db")
        if path.exists("tests/reference_annotation.db"):
            remove("tests/reference_annotation.db")
        if path.exists("tests/query_annotation1.db"):
            remove("tests/query_annotation1.db")
        if path.exists("tests/query_annotation2.db"):
            remove("tests/query_annotation2.db")
