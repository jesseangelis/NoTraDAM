import sqlite3
from os import path, remove

import numpy as np
import pytest

from src.analysis import Analysis
from src.database_handler import AnalysisDB, AnnotationDB


class TestAnalysis:
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
    def generate_transcript_ids(self):
        yield {
            "tests/reference_annotation.gtf": [
                "ENST111",
                "ENST122",
                "ENST212",
                "ENST221",
                "ENST231",
            ],
            "tests/query_annotation1.gtf": [
                "ENST121",
                "ENST122",
                "ENST212",
                "ENST231",
                "ENSTWrong",
            ],
            "tests/query_annotation2.gtf": [
                "ENST111",
                "ENST122",
                "ENST211",
                "ENST212",
                "ENST221",
                "ENST231",
            ],
        }

    @pytest.fixture
    def reference_db_path(self, generte_reference_gtf, generate_transcript_ids, dummy_args_instance):
        if path.exists("tests/reference_annotation.db"):
            remove("tests/reference_annotation.db")
        db_path = AnnotationDB.create_annotation_db(
            "tests/reference_annotation.gtf",
            "tests",
            generate_transcript_ids["tests/reference_annotation.gtf"],
            dummy_args_instance
        )
        yield db_path
        if path.exists(db_path):
            remove(db_path)

    @pytest.fixture
    def query_connection1_path(self, generte_query1_gtf, generate_transcript_ids, dummy_args_instance):
        if path.exists("tests/query_annotation1.db"):
            remove("tests/query_annotation1.db")
        db_path = AnnotationDB.create_annotation_db(
            "tests/query_annotation1.gtf",
            "tests",
            generate_transcript_ids["tests/query_annotation1.gtf"],
            dummy_args_instance
        )
        yield db_path
        if path.exists(db_path):
            remove(db_path)

    @pytest.fixture
    def query_connection2_path(self, generte_query2_gtf, generate_transcript_ids, dummy_args_instance):
        if path.exists("tests/query_annotation2.db"):
            remove("tests/query_annotation2.db")
        db_path = AnnotationDB.create_annotation_db(
            "tests/query_annotation2.gtf",
            "tests",
            generate_transcript_ids["tests/query_annotation2.gtf"],
            dummy_args_instance
        )
        yield db_path
        if path.exists(db_path):
            remove(db_path)

    @pytest.fixture
    def analysis_db_path(
        self, reference_db_path, query_connection1_path, query_connection2_path
    ):
        if path.exists("tests/analysis.db"):
            remove("tests/analysis.db")
        db_path = AnalysisDB.create_analysis_db(
            "tests", reference_db_path, [query_connection1_path, query_connection2_path]
        )
        yield db_path
        if path.exists(db_path):
            remove(db_path)

    @pytest.fixture
    def chromosome_sequences(self):
        yield {
            "chr1": "ATCTCAGGTGTACCGTTTCCGCCCCAACGGGTAGAAAGCCTGCAATCCCGCATCATTGATTGATGAAGTTCTATTAACTAGCCTAAACTACTGGGTACCGCTACAACTGTTAGACGCCGTTGTTATGATTCTGTCTTTCTGACATGAGAAATCCAGTGGGTAAGAGTGCTCGAGTAGCTTCATGTTGCTTTGGACTATCATCGATTCACGAGGGTGGTTGGCTTTTGCAACGCATCTCTCTAGCAAAATAAGTCGAGACTGCGCAGTGCGCCCCGTAAGTATTTCACGGGGAAGAAATACTTATGCGGGACCCGCTCACATCCTTCGGCTTGTAACTTACGAGTATACATAAAGTCTACGCAGACAGAGTTTAGGAACACACCGATGGTGACGGTGAG",
            "chr2": "AACCTCCTTCGGATTTGCGTCTGTGCCAATCAAGGCTTCTTCGCATGATCACTAAATGCAGATCCTTGTTGCGTATAGCACAAGGACATGAAGCCCTGCTAGGTAGGATGTAGCACTTGAACGCAGAGTTTGACAAAATGGTCGGGTCGAAGCGTACCAACAGATTTTATCGCGACCTGCTCGTGAGTCAATCCCTATACGGGAAGACGCAGCGAGATAAGGCCTCGCGTGTGCGAGCCCAGTTATATGAAGCACCTGTGCTCCAACTCAATCGCCGGTTGTGGGTTATGATTCCTCTCATCTTCTCACAGCGACTGAACAGCAGTTAGAATCAGGAGTAGAGGCGCCCGGATGTCCTTCCCGCGATAACACCGACATATCCTTACCATCTCCGGAGGAGCACAAGGAATCAGTGGGAAATAGAGACGCGCGCCAGGTAAGACGGTTCGGCAGAAACTATTCACTTCCAGAAGACACAACTACGAGCGTAGCGTACTGT",
        }

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

    @pytest.fixture
    def filled_analysis_db_path(
        self,
        reference_db_path,
        query_connection1_path,
        query_connection2_path,
        analysis_db_path,
        chromosome_sequences,
        dummy_args_instance,
    ):
        Analysis.compare_annotations(
            reference_db_path,
            query_connection1_path,
            analysis_db_path,
            chromosome_sequences,
            dummy_args_instance,
        )
        Analysis.compare_annotations(
            reference_db_path,
            query_connection2_path,
            analysis_db_path,
            chromosome_sequences,
            dummy_args_instance,
        )
        yield analysis_db_path
        if path.exists(analysis_db_path):
            remove(analysis_db_path)

    def test_compare_annotations(
        self,
        reference_db_path,
        query_connection1_path,
        query_connection2_path,
        analysis_db_path,
        chromosome_sequences,
        dummy_args_instance,
    ):
        """
        Description:
            Test the compare_annotations method of Analysis class.
            Compares the annotations of two query databases with a reference database.
            Asserts that the results match the expected values.

        Args:
            reference_db_path (str):
                path to the reference database.
            query_connection1_path (str):
                path to the first query database.
            query_connection2_path (str):
                path to the second query database.
            analysis_db_path (str):
                path to the analysis database.
            chromosome_sequences (dict):
                dictionary containing chromosome sequences.
            dummy_args_instance (argparse.ArgumentParser instance):
                ArgumentParser instance containing the arguments for the analysis.
        """
        # Compare annotations of the two query databases with the reference database.
        Analysis.compare_annotations(
            reference_db_path,
            query_connection1_path,
            analysis_db_path,
            chromosome_sequences,
            dummy_args_instance,
        )
        Analysis.compare_annotations(
            reference_db_path,
            query_connection2_path,
            analysis_db_path,
            chromosome_sequences,
            dummy_args_instance,
        )

        # Connect to the analysis database.
        analysis_db_connection = sqlite3.connect(analysis_db_path)

        # Assert that the results match the expected values.
        result1 = analysis_db_connection.execute(
            "SELECT transcript_id, query_annotation1, query_annotation2 FROM max_similarity"
        ).fetchall()
        expected1 = [
            ("ENST111", 0.0, 1.0),
            ("ENST122", 1.0, 1.0),
            ("ENST212", 1.0, 1.0),
            ("ENST221", 0.0, 1.0),
            ("ENST231", 1.0, 1.0),
        ]
        assert result1 == expected1
        result2 = analysis_db_connection.execute(
            "SELECT * FROM query_annotation1_transcripts"
        ).fetchall()
        expected2 = [
            ("ENST121", 4 / 5),
            ("ENST122", 1.0),
            ("ENST212", 1.0),
            ("ENST231", 1.0),
            ("ENSTWrong", 0.0),
        ]
        assert result2 == expected2
        result3 = analysis_db_connection.execute(
            "SELECT * FROM query_annotation2_transcripts"
        ).fetchall()
        expected3 = [
            ("ENST111", 1.0),
            ("ENST122", 1.0),
            ("ENST211", 3 / 5),
            ("ENST212", 1.0),
            ("ENST221", 1.0),
            ("ENST231", 1.0),
        ]
        assert result3 == expected3

        # Close the connection to the analysis database.
        analysis_db_connection.close()

    def test_calculate_metrics(self, filled_analysis_db_path):
        """
        Description:
            Test the calculate_metrics method of the Analysis class.

        Args:
            filled_analysis_db_path (str): 
                The path to the filled analysis database.

        Returns:
            None
        """
        # Define the query names to be used in the test
        query_names = ["query_annotation1", "query_annotation2"]

        # Call the calculate_metrics method with the filled analysis database and query names
        Analysis.calculate_metrics(filled_analysis_db_path, query_names)

        # Connect to the filled analysis database
        filled_analysis_db = sqlite3.connect(filled_analysis_db_path)

        # Retrieve the results for query_annotation1 and query_annotation2 from the filled analysis database
        result_query_annotation1 = filled_analysis_db.execute(
            "SELECT * FROM query_annotation1_results"
        ).fetchall()
        result_query_annotation2 = filled_analysis_db.execute(
            "SELECT * FROM query_annotation2_results"
        ).fetchall()

        # Define the expected results for query_annotation1 and query_annotation2
        expected_query_annotation1 = list(
            zip(
                list(np.linspace(0.81, 1.0, num=20))[::-1],
                [3] * 20,
                [3] * 20,
                [2] * 20,
                [2] * 20,
                [2 / 5] * 20,
                [3 / 5] * 20,
            )
        ) + list(
            zip(
                list(np.linspace(0.75, 0.8, num=6))[::-1],
                [4] * 6,
                [3] * 6,
                [1] * 6,
                [2] * 6,
                [1 / 5] * 6,
                [3 / 5] * 6,
            )
        )
        expected_query_annotation2 = list(
            zip(
                list(np.linspace(0.75, 1, num=26))[::-1],
                [5] * 26,
                [5] * 26,
                [1] * 26,
                [0] * 26,
                [1 / 6] * 26,
                [1] * 26,
            )
        )

        # Check that the actual results match the expected results for query_annotation1 and query_annotation2
        for r, e in zip(result_query_annotation1, expected_query_annotation1):
            assert np.allclose(r, e, atol=0.00001), (r, e)
        for r, e in zip(result_query_annotation2, expected_query_annotation2):
            assert np.allclose(r, e, atol=0.00001), (r, e)

        # Close the filled analysis database connection
        filled_analysis_db.close()
