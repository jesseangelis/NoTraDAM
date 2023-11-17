import random
from os import path, remove

import pytest

from src.database_handler import AnnotationDB
from src.sequence_handler import FastaParser, Sequences


class TestFastaParser:
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

    def test_parse_fasta(self, generate_fasta_file):
        fasta_path = "tests/fasta.fasta"
        expected_output = {
            "chr1": "ATCTCAGGTGTACCGTTTCCGCCCCAACGGGTAGAAAGCCTGCAATCCCGCATCATTGATTGATGAAGTTCTATTAACTAGCCTAAACTACTGGGTACCGCTACAACTGTTAGACGCCGTTGTTATGATTCTGTCTTTCTGACATGAGAAATCCAGTGGGTAAGAGTGCTCGAGTAGCTTCATGTTGCTTTGGACTATCATCGATTCACGAGGGTGGTTGGCTTTTGCAACGCATCTCTCTAGCAAAATAAGTCGAGACTGCGCAGTGCGCCCCGTAAGTATTTCACGGGGAAGAAATACTTATGCGGGACCCGCTCACATCCTTCGGCTTGTAACTTACGAGTATACATAAAGTCTACGCAGACAGAGTTTAGGAACACACCGATGGTGACGGTGAG",
            "chr2": "AACCTCCTTCGGATTTGCGTCTGTGCCAATCAAGGCTTCTTCGCATGATCACTAAATGCAGATCCTTGTTGCGTATAGCACAAGGACATGAAGCCCTGCTAGGTAGGATGTAGCACTTGAACGCAGAGTTTGACAAAATGGTCGGGTCGAAGCGTACCAACAGATTTTATCGCGACCTGCTCGTGAGTCAATCCCTATACGGGAAGACGCAGCGAGATAAGGCCTCGCGTGTGCGAGCCCAGTTATATGAAGCACCTGTGCTCCAACTCAATCGCCGGTTGTGGGTTATGATTCCTCTCATCTTCTCACAGCGACTGAACAGCAGTTAGAATCAGGAGTAGAGGCGCCCGGATGTCCTTCCCGCGATAACACCGACATATCCTTACCATCTCCGGAGGAGCACAAGGAATCAGTGGGAAATAGAGACGCGCGCCAGGTAAGACGGTTCGGCAGAAACTATTCACTTCCAGAAGACACAACTACGAGCGTAGCGTACTGT",
        }
        assert FastaParser.parse_fasta(fasta_path) == expected_output


class TestSequences:
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

    @pytest.fixture
    def gtf_file(self):
        content = """##description: reference In reduced not included ENST122, ENST212, ENST231
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
            f.write(content)
        yield "tests/reference_annotation.gtf"
        if path.exists("tests/reference_annotation.gtf"):
            remove("tests/reference_annotation.gtf")

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

    def test_get_sequence(self, gtf_file, chromosome_sequences, generate_fasta_file, dummy_args_instance):
        """
        Description:
            Test the get_sequence method of the Sequences class.

        Args:
            gtf_file (str): The path to the GTF file.
            chromosome_sequences (dict): A dictionary containing the chromosome sequences.
            generate_fasta_file (function): A function that generates a FASTA file.

        Returns:
            None
        """
        # Remove the reference annotation database if it exists
        if path.exists("tests/reference_annotation.db"):
            remove("tests/reference_annotation.db")

        # Create the annotation database
        db_path = AnnotationDB.create_annotation_db(gtf_file, "tests", None, dummy_args_instance)

        # Parse the FASTA file
        chromosome_sequences = FastaParser.parse_fasta("tests/fasta.fasta")

        # Test the get_sequence method with transcript_id = "ENST111"
        transcript_id = "ENST111"
        assert (
            Sequences.get_sequence(transcript_id, db_path, chromosome_sequences)
            == "ATC"
        )

        # Test the get_sequence method with transcript_id = "ENST211"
        transcript_id = "ENST211"
        assert (
            Sequences.get_sequence(transcript_id, db_path, chromosome_sequences)
            == "TACTA"
        )

        # Remove the reference annotation database if it exists
        if path.exists("tests/reference_annotation.db"):
            remove("tests/reference_annotation.db")

    def test_reverse_complement(self):
        """
        Description:
            Test the reverse complement function of the Sequences class.
        """
        # Test with a simple DNA sequence
        assert Sequences.reverse_complement("ATCG") == "CGAT"

        # Test with a more complex DNA sequence
        assert Sequences.reverse_complement("ATCGAT") == "ATCGAT"

    def test_compare_sequences(self):
        """
        Description:
            Test the compare_sequences method of the Sequences class.
            This test checks if the compare_sequences method returns the expected similarity score
            for two given DNA sequences and a set of alignment arguments.
        """
        # Define the alignment arguments and two DNA sequences to compare
        alignment_args = (1, 0, 0, 0)
        seq_1 = "ATCGAT"
        seq_2 = "ATCGAT"

        # Check if the similarity score for identical sequences is 1.0
        assert Sequences.compare_sequences(seq_1, seq_2, alignment_args) == 1.0

        # Change one nucleotide in the second sequence and check if the similarity score is 5/7
        seq_1 = "ATCGAT"
        seq_2 = "ATCAGA"
        assert Sequences.compare_sequences(seq_1, seq_2, alignment_args) == 5 / 7
