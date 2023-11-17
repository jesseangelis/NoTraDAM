# NoTraDAM
No<sub>vel</sub> Tra<sub>nscript</sub> D<sub>etection</sub> A<sub>lgorithm</sub> M<sub>etrics</sub>

A tool to evaluate at tool's capabilities to identify novel transcripts based on GTF annotaion files. The evaluation is based on true positiv rate (TPR) and false discovery rate (FDR) depending on a sequence similarity threshold.

## How to use

Clone repo

```
cd target_directory
git clone https://git.rz.tu-bs.de/scibiome/phd-students/lrs_isoid.git
cd NoTraDam
```

Install requirments
with pip:
```
pip install -r requirements.txt
``` 
Test with:
```
pytest
```

## Method

Before utilizing NoTraDAM, you need to follow these steps:
- Have a genome annotation file (reference annotation) and chromosome sequence file (genome fasta) ready
- Use these to simulate reads
- Remove a subest of the transcripts and their exons from the gtf file
- Use a novel transcript detection algortihm with the simulated reads and the reduced annotation file.
- The algorithm should produce a new annotaiton file (query annotation) where the known and novel transcripts are included

Run command:
```
python notradam.py -r reference_annoation.gtf -n query_annotation1.gtf query_annotation2.gtf -i transcript_id.txt -f genome_fasta.fasta
```

| Argument                            | Description                                                                                                      |
| ------------------------------------| ---------------------------------------------------------------------------------------------------------------- |
| `-r, --reference_annotation`         | The reference annotation as .gtf file                                                                            |
| `-n, --query_annotations`            | The predicted annotation or annotations as .gtf file/s                                                           |
| `-i, --transcript_ids`               | Text file, containing the transcript IDs, to be included in the analysis. If not given, all transcripts will be included |
| `-f, --genome_fasta`                 | Fasta file with chromosome sequences. The chromosome IDs must be the same as in the annotations               |
| `-o, --output_directory`             | Directory used for the outputs                                                                                  |
| `-k, --keep_temporary`               | If given, temporary files are kept. This may increase performance for reruns but may result in wrong results/error if files have the same paths but are not the same |
| `-m, --match_score`                  | Match score. Default = 1. See Biopython Align.PairwiseAlignment documentation                                   |
| `-x, --mismatch_score`               | Mismatch score. Default = 0. See Biopython Align.PairwiseAlignment documentation                                 |
| `-g, --open_gap_score`               | Open gap score. Default = 0. See Biopython Align.PairwiseAlignment documentation                                 |
| `-e, --extend_gap_score`             | Extended gap score. Default = 0. See Biopython Align.PairwiseAlignment documentation                              |
| `-t, --target_end_gap_score`         | Target end gap score. Default = 0. See Biopython Align.PairwiseAlignment documentation                            |
| `-q, --query_end_gap_score`          | Query end gap score. Default = 0. See Biopython Align.PairwiseAlignment documentation                              |
| `-v, --verbose`                      | If given, more information is printed to the console                                                             |
| `-c, --cores`                        | Number of cores used for the analysis. Default = 1                                                               |


**transcript_id** should be a plain text file. If given, only these transcripts are used for the analysis. This can be used to select only the reference transcript that were removed prior to running the detection tool and the query transcripts marked as novel. The structure should be as follows:

```
@path_to_reference_annotation.gtf
reference_transcript_id_1
reference_transcript_id_2 

@path_to_query_annotation.gtf
novel_transcript_id_1
novel_transcript_id_2 

@path_to_query_annotation.gtf
novel_transcript_id_1
novel_transcript_id_2

and so on...
```
NoTraDam works as follows:
- 1. For every chromosome-strand combination a transcript subset is selected from the reference and query annotation
- 2. Matches where there is an overlap of the transcript start and end postions are found
- 3. For each match the best allignment is found and a simularity score is calculated. This is the number of identical bases divided by the length of the allignment.
- 4. For each reference and query transcript the maximum found similarity is saved.
- 5. Based on a similarity threshold tpr and fdr are calculated


### True positive rate (Sensetivity)
The sensitivity measures the tool's ability to completely predict transcipts that is not in the annotaiton. It measures false negatives.

TPR $= \frac{TP}{TP + FN}$

Where TP are the number of reference transcripts that are not known to the query that have a sufficiant similiarity to a transcript with the query (larger or equal to the set threshold) and FN are the reference transcripts that are not known to the query that do not have a corresponding predicted transcript with sufficiant similiarity (lower then the threshold).


### False Discovery Rate 
The false discovery rate measures a tools' ability to only predict acutal novel transcripts. It measures false positives.

FDR $= \frac{FP}{TP + FP}$

Where TP are the number of query transcripts that have a corresponding reference transcript that was removed for the detection with sufficiant similiarity (larger or equal to the set threshold) and FP are the number of query transcripts that no not have a corresponding removed transcript with sufficiant similiarity (lower then the threshold).

# Output

The output directory is structured as follows:

```
output_dir/
│
└───tempdir/
│   │   reference_annotation.db
│   │   query1_annotation.db
│   │   query2_annotation.db
│   │   ...
│
│   analysis.db
│   result_combined.csv
│   query1_result.csv
│   query2_result.csv
│   ...
```

<div style="background-color: #000; color: #fff; border: 1px solid #333; padding: 10px; margin-bottom: 0;">
    <strong>Note:</strong>
    <p style="margin-bottom: 0;">NoTraDAM is not suitable for trans-splicing events as the exons for a transcript are ordered by their start position.</p>
</div>


creator: Jesse Angelis (jesse@angelis.info)