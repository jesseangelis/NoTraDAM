import sqlite3
from time import process_time

from Bio import Align, SeqIO


class FastaParser:
    def parse_fasta(path: str) -> dict():
        """
        Parses the chromosome fasta file and creates a dict with chromosome IDs as keys and sequences as values

        Parameters
        -----------
        path : str
            path to chromosme fasta file

        Returns
        -----------
        chromosome_sequences : dict
            Dictionary with chromosome names as keys and seqence strings as values
        """
        chromosome_sequences = {}
        records = SeqIO.parse(open(path), "fasta")
        for chromosome in records:
            chromosome_sequences[chromosome.id] = chromosome.seq
        return chromosome_sequences


class Sequences:
    def get_sequence(
        transcript_id: str, db_path: str, chromosome_sequences: dict) -> str:
        """
        Get the sequence of a transcript.

        Parameters
        ----------
        transcript_id : str
            The ID of the transcript.
        connection : sqlite3.Connection
            A connection to the annotation database.
        chromosome_sequences : dict of str : str
            A dictionary mapping chromosome names to their sequences.

        Returns
        -------
        str
            The sequence of the transcript.
        """
        start_time = process_time()
        while True:
            try:
                connection = sqlite3.connect(db_path)
                strand = connection.execute(f"SELECT strand FROM exons WHERE transcript_id = ?", (transcript_id,)).fetchone()[0]
                exons = connection.execute(
                    f"SELECT * FROM exons WHERE transcript_id = ? ORDER BY CASE WHEN ? = '+' THEN start END ASC, CASE WHEN ? = '-' THEN start END DESC",
                    (transcript_id, strand, strand),
                ).fetchall()
                chromosome = exons[0][1]
                sequence = ""
                for exon in exons:
                    if strand == "+":
                        start = exon[3]
                        end = exon[4]
                        sequence += chromosome_sequences[chromosome][start - 1 : end]
                    elif strand == "-":
                        start = exon[4]
                        end = exon[3]
                        sequence += Sequences.reverse_complement(
                            chromosome_sequences[chromosome][start - 1 : end]
                        )
                connection.close()
            except sqlite3.Error:
                if process_time() - start_time > 3:
                    raise sqlite3.OperationalError("Timeout due to database lock")
            else:
                break
        return sequence

    def reverse_complement(sequence: str) -> str:
        """
        Calculates the reverse complement of a sequence.

        Parameters
        ----------
        sequence : str
            The sequence to be reversed.

        Returns
        -------
        str
            The reverse complement of the sequence.
        """
        complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
        return "".join(complement[base] for base in reversed(sequence))

    def compare_sequences(
        reference_sequence: str, query_sequence: str, alignment_args
    ) -> float:
        """
        Compares two DNA sequences and returns the best alignment.

        Parameters
        ----------
        reference_sequence : str
            The reference DNA sequence.
        query_sequence : str
            The query DNA sequence.
        args : argparse.ArgumentParser object
            The parsed command-line arguments.

        Returns
        -------
        float
            The similarity score of the best alignment.
        """
        aligner = Align.PairwiseAligner()
        (
            aligner.match_score,
            aligner.mismatch_score,
            aligner.open_gap_score,
            aligner.extend_gap_score,
        ) = alignment_args
        alignment = aligner.align(reference_sequence, query_sequence)[0]
        count = 0
        for s1, s2 in zip(alignment[0], alignment[1]):
            if s1 == s2:
                count += 1
        similarity = count / len(alignment[0])
        return similarity
