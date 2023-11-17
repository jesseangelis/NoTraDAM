from os import remove
from os.path import exists
from shutil import rmtree

import numpy as np
import pytest
from Bio.Seq import Seq

from src.main import Pipeline


class TestPipeline:
    @pytest.fixture(autouse=True)
    def dummy_pipeline_instance(self):
        class DummyPipeline(Pipeline):
            def __init__(self):
                pass

        yield DummyPipeline()

    @pytest.fixture(autouse=True)
    def dummy_args_instance(self):
        class DummyArgs:
            def __init__(self):
                self.reference_annotation = "tests/reference_annotation.gtf"
                self.query_annotations = [
                    "tests/query_annotation1.gtf",
                    "tests/query_annotation2.gtf",
                ]
                self.transcript_ids = "tests/transcript_ids.txt"
                self.genome_fasta = "tests/fasta.fasta"
                self.output_directory = "tests/output"
                self.keep_temporary = True
                self.verbose = True
                self.cores = 1
                self.match_score = 1
                self.mismatch_score = 0
                self.open_gap_score = 0
                self.extend_gap_score = 0

        yield DummyArgs()

    @pytest.fixture
    def transcript_ids_file(self):
        file_path = "tests/transcript_ids.txt"
        with open(file_path, "w") as f:
            f.write("@tests/reference_annotation.gtf\n")
            f.write("ENST111\nENST122\nENST211\nENST212\n ENST221\nENST231\n")
            f.write("@tests/query_annotation1.gtf\n")
            f.write("ENST111\nENST122\nENST211\nENST212\nENST221\n")
            f.write("\n")
            f.write("@tests/query_annotation2.gtf\n")
            f.write("ENST111\nENST122\nENST211\nENST212\nENST231\n")
        yield file_path
        remove(file_path)

    @pytest.fixture
    def generate_gtf_files(self):
        reference_annotation_content = """##description: reference In reduced not included ENST111, ENST122, ENST211, ENST212, ENST221, ENST231
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
        query_annotation1_content = """##description: query 1 Found: ENST122 False: ENST212 Not found: ENST231
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
chr2	TEST	exon	103	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST211"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 2; exon_id "Exon2112";
chr2	TEST	transcript	103	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript212";
chr2	TEST	exon	103	100	.	-	.	gene_id "ENSG21"; transcript_id "ENST212"; gene_type "lncRNA"; gene_name "Gene21"; transcript_type "lncRNA"; transcript_name "Transcript211"; exon_number 1; exon_id "Exon2112";
chr2	TEST	gene	300	302	.	+	.	gene_id "ENSG22"; gene_type "lncRNA"; gene_name "Gene22";
chr2	TEST	transcript	300	302	.	+	.	gene_id "ENSG22"; transcript_id "ENST221"; gene_type "lncRNA"; gene_name "Gene22"; transcript_type "lncRNA"; transcript_name "Transcript221";
chr2	TEST	exon	300	302	.	+	.	gene_id "ENSG22"; transcript_id "ENST221"; gene_type "lncRNA"; gene_name "Gene22"; transcript_type "lncRNA"; transcript_name "Transcript221"; exon_number 1; exon_id "Exon2211";"""
        query_annotation2_content = """##description: query 2 Not found: ENST122, Found: ENST212, ENST231
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
        with open("tests/query_annotation1.gtf", "w") as f:
            f.write(query_annotation1_content)
        with open("tests/query_annotation2.gtf", "w") as f:
            f.write(query_annotation2_content)
        yield
        if exists("tests/reference_annotation.gtf"):
            remove("tests/reference_annotation.gtf")
        if exists("tests/query_annotation1.gtf"):
            remove("tests/query_annotation1.gtf")
        if exists("tests/query_annotation2.gtf"):
            remove("tests/query_annotation2.gtf")

    @pytest.fixture
    def generate_fasta_file(self):
        with open("tests/fasta.fasta", "w") as f:
            f.write(">chr1\n")
            f.write(
                "ATCTCAGGTGTACCGTTTCCGCCCCAACGGGTAGAAAGCCTGCAATCCCGCATCATTGATTGATGAAGTTCTATTAACTAGCCTAAACTACTGGGTACCGCTACAACTGTTAGACGCCGTTGTTATGATTCTGTCTTTCTGACATGAGAAATCCAGTGGGTAAGAGTGCTCGAGTAGCTTCATGTTGCTTTGGACTATCATCGATTCACGAGGGTGGTTGGCTTTTGCAACGCATCTCTCTAGCAAAATAAGTCGAGACTGCGCAGTGCGCCCCGTAAGTATTTCACGGGGAAGAAATACTTATGCGGGACCCGCTCACATCCTTCGGCTTGTAACTTACGAGTATACATAAAGTCTACGCAGACAGAGTTTAGGAACACACCGATGGTGACGGTGAG\n"
            )
            f.write(">chr2\n")
            f.write(
                "AACCTCCTTCGGATTTGCGTCTGTGCCAATCAAGGCTTCTTCGCATGATCACTAAATGCAGATCCTTGTTGCGTATAGCACAAGGACATGAAGCCCTGCTAGGTAGGATGTAGCACTTGAACGCAGAGTTTGACAAAATGGTCGGGTCGAAGCGTACCAACAGATTTTATCGCGACCTGCTCGTGAGTCAATCCCTATACGGGAAGACGCAGCGAGATAAGGCCTCGCGTGTGCGAGCCCAGTTATATGAAGCACCTGTGCTCCAACTCAATCGCCGGTTGTGGGTTATGATTCCTCTCATCTTCTCACAGCGACTGAACAGCAGTTAGAATCAGGAGTAGAGGCGCCCGGATGTCCTTCCCGCGATAACACCGACATATCCTTACCATCTCCGGAGGAGCACAAGGAATCAGTGGGAAATAGAGACGCGCGCCAGGTAAGACGGTTCGGCAGAAACTATTCACTTCCAGAAGACACAACTACGAGCGTAGCGTACTGT\n"
            )
        yield
        remove("tests/fasta.fasta")

    def test_create_output_dir(self, dummy_pipeline_instance, dummy_args_instance):
        """
        Description:
            Test case to check if the output directory is created successfully.

        Args:
            dummy_pipeline_instance: An instance of the pipeline class.
            dummy_args_instance: An instance of the arguments class.

        Returns:
            None
        """
        # Remove the output directory if it already exists
        if exists(dummy_args_instance.output_directory):
            rmtree(dummy_args_instance.output_directory)

        # Create the output directory
        dummy_pipeline_instance.create_output_dir(dummy_args_instance)

        # Check if the output directory and temp directory are created
        assert exists(dummy_args_instance.output_directory)
        assert exists(dummy_args_instance.output_directory + "/temp")

        # Remove the output directory after the test
        if exists(dummy_args_instance.output_directory):
            rmtree(dummy_args_instance.output_directory)

    def test_parse_transcript_ids(
        self, dummy_pipeline_instance, dummy_args_instance, transcript_ids_file
    ):
        """
        Description:
            Test the parse_transcript_ids method of the Pipeline class.
            This method should parse the transcript IDs from the input GTF files and store them in the
            transcript_ids attribute of the Pipeline instance.

        Args:
            dummy_pipeline_instance (Pipeline): An instance of the Pipeline class for testing.
            dummy_args_instance (Namespace): An instance of the argparse.Namespace class for testing.
            transcript_ids_file (str): The path to the transcript IDs file for testing.

        Returns:
            None
        """
        # Call the parse_transcript_ids method of the Pipeline instance
        dummy_pipeline_instance.parse_transcript_ids(dummy_args_instance)

        # Define the expected transcript IDs for each GTF file
        expected = {
            "tests/reference_annotation.gtf": np.array(
                ["ENST111", "ENST122", "ENST211", "ENST212", "ENST221", "ENST231"]
            ),
            "tests/query_annotation1.gtf": np.array(
                ["ENST111", "ENST122", "ENST211", "ENST212", "ENST221"]
            ),
            "tests/query_annotation2.gtf": np.array(
                ["ENST111", "ENST122", "ENST211", "ENST212", "ENST231"]
            ),
        }

        # Get the actual transcript IDs from the Pipeline instance
        actual = dummy_pipeline_instance.transcript_ids

        # Check that the actual transcript IDs match the expected transcript IDs for each GTF file
        assert np.array_equal(
            actual["tests/reference_annotation.gtf"],
            expected["tests/reference_annotation.gtf"],
        )
        assert np.array_equal(
            actual["tests/query_annotation1.gtf"],
            expected["tests/query_annotation1.gtf"],
        )
        assert np.array_equal(
            actual["tests/query_annotation2.gtf"],
            expected["tests/query_annotation2.gtf"],
        )

    def test_create_annotation_dbs(
        self,
        dummy_pipeline_instance,
        dummy_args_instance,
        generate_gtf_files,
        transcript_ids_file,
    ):
        """
        Description:
            Test case to check if the annotation databases are created successfully.

        Args:
            dummy_pipeline_instance: An instance of the pipeline class.
            dummy_args_instance: An instance of the argparse class.
            generate_gtf_files: A fixture to generate GTF files.
            transcript_ids_file: A fixture to generate transcript IDs file.
        """
        # Remove existing annotation databases if any
        if exists(dummy_args_instance.output_directory + "/temp/query_annotation1.db"):
            remove(dummy_args_instance.output_directory + "/temp/query_annotation1.db")
        if exists(dummy_args_instance.output_directory + "/temp/query_annotation2.db"):
            remove(dummy_args_instance.output_directory + "/temp/query_annotation2.db")
        if exists(
            dummy_args_instance.output_directory + "/temp/reference_annotation.db"
        ):
            remove(
                dummy_args_instance.output_directory + "/temp/reference_annotation.db"
            )

        # Create output directory and parse transcript IDs
        dummy_pipeline_instance.create_output_dir(dummy_args_instance)
        dummy_pipeline_instance.parse_transcript_ids(dummy_args_instance)

        # Create annotation databases
        dummy_pipeline_instance.create_annotation_dbs(dummy_args_instance)

        # Assertions to check if the databases are created successfully
        assert isinstance(dummy_pipeline_instance.reference_annotation_db_path, str)
        assert isinstance(
            dummy_pipeline_instance.query_annotation_db_paths[
                dummy_args_instance.query_annotations[0]
            ],
            str,
        )
        assert isinstance(
            dummy_pipeline_instance.query_annotation_db_paths[
                dummy_args_instance.query_annotations[1]
            ],
            str,
        )
        assert exists(
            dummy_args_instance.output_directory + "/temp/query_annotation1.db"
        )
        remove(dummy_args_instance.output_directory + "/temp/query_annotation1.db")
        assert exists(
            dummy_args_instance.output_directory + "/temp/query_annotation2.db"
        )
        remove(dummy_args_instance.output_directory + "/temp/query_annotation2.db")
        assert exists(
            dummy_args_instance.output_directory + "/temp/reference_annotation.db"
        )
        remove(dummy_args_instance.output_directory + "/temp/reference_annotation.db")

    def test_get_chromosome_sequences(
        self, dummy_pipeline_instance, dummy_args_instance, generate_fasta_file
    ):
        """
        Description:
            Test the get_chromosome_sequences method of Pipeline class.
            This method should retrieve the chromosome sequences from the provided fasta file
            and store them in a dictionary attribute of the Pipeline instance.

        Args:
            dummy_pipeline_instance (Pipeline): An instance of the Pipeline class.
            dummy_args_instance (Namespace): An instance of the argparse.Namespace class.
            generate_fasta_file (fixture): A fixture that generates a temporary fasta file.

        Returns:
            None

        Raises:
            AssertionError: If any of the test conditions are not met.
        """
        # Call the method being tested
        dummy_pipeline_instance.get_chromosome_sequences(dummy_args_instance)

        # Check that the chromosome_sequences attribute is a dictionary
        assert type(dummy_pipeline_instance.chromosome_sequences) == dict

        # Check that the chromosome_sequences attribute has the expected number of entries
        assert len(dummy_pipeline_instance.chromosome_sequences) == 2

        # Check that the chromosome sequences are of the expected type
        assert isinstance(dummy_pipeline_instance.chromosome_sequences["chr1"], Seq)
        assert isinstance(dummy_pipeline_instance.chromosome_sequences["chr2"], Seq)

    def test_create_analysis_db(
        self,
        dummy_pipeline_instance,
        dummy_args_instance,
        generate_gtf_files,
        transcript_ids_file,
    ):
        """
        Description:
            Test case to check if the analysis database is created successfully.

        Args:
            dummy_pipeline_instance: An instance of the pipeline class.
            dummy_args_instance: An instance of the arguments class.
            generate_gtf_files: A fixture to generate GTF files.
            transcript_ids_file: A fixture to generate transcript IDs file.

        Returns:
            None
        """
        if exists(dummy_args_instance.output_directory + "/analysis.db"):
            remove(dummy_args_instance.output_directory + "/analysis.db")
        dummy_pipeline_instance.create_output_dir(dummy_args_instance)
        dummy_pipeline_instance.parse_transcript_ids(dummy_args_instance)
        dummy_pipeline_instance.create_annotation_dbs(dummy_args_instance)
        dummy_pipeline_instance.create_analysis_db(dummy_args_instance)
        assert isinstance(dummy_pipeline_instance.analysis_db_path, str)
        assert exists(dummy_args_instance.output_directory + "/analysis.db")
        remove(dummy_args_instance.output_directory + "/analysis.db")
        rmtree(dummy_args_instance.output_directory)

    def test_compare_reference_with_queries(
        self,
        dummy_pipeline_instance,
        dummy_args_instance,
        generate_gtf_files,
        generate_fasta_file,
        transcript_ids_file,
    ):
        """
        Description:
            Test the compare_reference_with_queries method of the Pipeline class.
            This method compares the reference genome with the queries and generates an analysis database.
            The method also checks if the analysis database path is a string and removes the analysis database and output directory.

        Args:
            dummy_pipeline_instance (Pipeline): An instance of the Pipeline class.
            dummy_args_instance (Namespace): An instance of the Namespace class.
            generate_gtf_files (fixture): A fixture that generates GTF files.
            generate_fasta_file (fixture): A fixture that generates a FASTA file.
            transcript_ids_file (fixture): A fixture that generates a transcript IDs file.
        """
        # Remove analysis database if it already exists
        if exists(dummy_args_instance.output_directory + "/analysis.db"):
            remove(dummy_args_instance.output_directory + "/analysis.db")

        # Get chromosome sequences
        dummy_pipeline_instance.get_chromosome_sequences(dummy_args_instance)

        # Create output directory
        dummy_pipeline_instance.create_output_dir(dummy_args_instance)

        # Parse transcript IDs
        dummy_pipeline_instance.parse_transcript_ids(dummy_args_instance)

        # Create annotation databases
        dummy_pipeline_instance.create_annotation_dbs(dummy_args_instance)

        # Create analysis database
        dummy_pipeline_instance.create_analysis_db(dummy_args_instance)

        # Compare reference with queries and generate analysis database
        dummy_pipeline_instance.compare_reference_with_queries(dummy_args_instance)

        # Check if analysis database path is a string
        assert isinstance(dummy_pipeline_instance.analysis_db_path, str)

        # Remove analysis database
        remove(dummy_pipeline_instance.analysis_db_path)

        # Remove output directory
        rmtree(dummy_args_instance.output_directory)

    def test_run(
        self,
        dummy_args_instance,
        generate_gtf_files,
        generate_fasta_file,
        transcript_ids_file,
    ):
        """
        Description:
            Test the run method of the Pipeline class.

        Args:
            dummy_args_instance (argparse.Namespace): An instance of argparse.Namespace
                containing the command line arguments.
            generate_gtf_files (function): A fixture that generates GTF files.
            generate_fasta_file (function): A fixture that generates a FASTA file.
            transcript_ids_file (function): A fixture that generates a transcript IDs file.
        """
        # Remove analysis.db if it already exists
        if exists(dummy_args_instance.output_directory + "/analysis.db"):
            remove(dummy_args_instance.output_directory + "/analysis.db")

        # Create a Pipeline instance and check if analysis_db_path is a string
        p = Pipeline(dummy_args_instance)
        assert isinstance(p.analysis_db_path, str)

        # Check if analysis.db was created
        assert exists(dummy_args_instance.output_directory + "/analysis.db")

        # Remove the output directory
        rmtree(dummy_args_instance.output_directory)
