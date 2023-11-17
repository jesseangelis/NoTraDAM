import shutil
import subprocess
from os import remove

import pytest


class TestNotradam:
    @pytest.fixture
    def subprocess_script(self):
        shutil.copyfile("notradam.py", "tests/subprocess_script.py")
        with open("tests/subprocess_script.py", "r") as f:
            lines = f.readlines()
        with open("tests/subprocess_script.py", "w") as f:
            for line in lines:
                if "from src.main import Pipeline" in line:
                    continue
                if 'if __name__ == "__main__":' in line:
                    break
                f.write(line)
            f.write("try:\n")
            f.write("\targs = get_arguments()\n")
            f.write(
                "\tprint(args.reference_annotation, args.query_annotations, args.transcript_ids, args.genome_fasta, args.output_directory, args.keep_temporary , args.verbose, args.cores)\n"
            )
            f.write("except Exception as e:\n")
            f.write("\tprint(e)")
        yield
        remove("tests/subprocess_script.py")

    def test_get_arguments(self, subprocess_script):
        """
        Description:
            Test the get_arguments function of the subprocess_script.py script.

        Args:
            subprocess_script: A fixture that returns the path to the subprocess_script.py script.

        Returns:
            None
        """
        # Define input arguments
        reference_annotation = "reference_annotation.gtf"
        query_annotation1 = "query_annotation1.gtf"
        query_annotation2 = "query_annotation2.gtf"
        transcript_ids = "transcript_ids.txt"
        genome_fasta = "fasta.fa"
        output_directory = "output_directory"
        keep_temporary = True
        verbose = True
        cores = "4"

        # Define subprocess command
        subprocess_command = [
            "python",
            "tests/subprocess_script.py",
            "-r",
            reference_annotation,
            "-n",
            query_annotation1,
            query_annotation2,
            "-i",
            transcript_ids,
            "-f",
            genome_fasta,
            "-o",
            output_directory,
            "-k",
            "-v",
            "-c",
            cores,
        ]

        # Run subprocess command and capture output
        output = subprocess.run(subprocess_command, capture_output=True, text=True)

        # Define expected output
        expected = f"{reference_annotation} ['{query_annotation1}', '{query_annotation2}'] {transcript_ids} {genome_fasta} {output_directory} {keep_temporary} {verbose} {int(cores)}"

        # Assert that actual output matches expected output
        assert output.stdout.strip() == expected, output.stdout.strip()
