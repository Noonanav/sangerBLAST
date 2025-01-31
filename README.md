# Sanger Sequencing Data and 16S Sequence Processor

This script processes Sanger sequencing data in a 2-step process: 

1. Paired forward and reverse reads are aligned using the `pairwise2.align.globalms` function to generate consensus sequences
2. A BLAST search against specified databases (either local and web-based) and top matches are reported in the output file.

This script can also handle existing 16S sequences in FASTA format using `--input_type fasta`

## Requirements

- Python 3.x
- Biopython
- pandas

## Installation

To install the required Python packages, you can use `pip` to install in a `conda` environment:

```
conda create -n sangerBLAST python=3
conda activate sangerBLAST

pip install biopython pandas
```

## Usage

To run the script, use the following command:

```
sangerBLAST.py [-h] [--output OUTPUT] [--consensus_fasta CONSENSUS_FASTA] [--trimmed_fasta TRIMMED_FASTA] [--verbose] [--web_blast] [--blast_Ns] [--output_dir OUTPUT_DIR] [--input_type {ab1,fasta}] input_directory

Process Sanger sequencing data, align reads, generate consensus sequences, and perform BLAST.

positional arguments:
  input_directory       Directory containing .ab1 files

optional arguments:
  -h, --help            show this help message and exit
  --output OUTPUT       Output file for BLAST results
  --consensus_fasta CONSENSUS_FASTA
                        Output file for consensus sequences in FASTA format
  --trimmed_fasta TRIMMED_FASTA
                        Output file for trimmed consensus sequences in FASTA format
  --verbose             Enable verbose output
  --web_blast           BLAST against NCBI database instead of local BLAST
  --blast_Ns            BLAST concensus sequences with Ns
  --output_dir OUTPUT_DIR
                        Directory to place output files
  --input_type {ab1,fasta}
                        Type of input data (default: ab1)
```

### Arguments

| Argument            | Description                                                                    | Default                           |
|---------------------|--------------------------------------------------------------------------------|-----------------------------------|
| `--output`          | Output file for BLAST results                                                  | `blast_results.csv`                 |
| `--output_dir`      | Directory to save output files                                                 | `.`                                 |
| `--consensus_fasta` | Output file for consensus sequences in FASTA format                            | `consensus_sequences.fasta`         |
| `--trimmed_fasta`   | Output file for trimmed consensus sequences in FASTA format                    | `trimmed_consensus_sequences.fasta` |
| `--verbose`         | Enable verbose output                                                          |         `False`                     |
| `--web_blast`       | BLAST against NCBI database instead of local BLAST                             |     `False`                         |
| `--blast_Ns`         | BLAST concensus sequences with Ns                                           |      `False`                        |
| `--input_type`      | Type of input data (options: ab1 / fasta)                                      | `ab1`                               |    

### Input

The script assumes the following input format:

| Input            | Description                                                                                      |
|-----------------------|--------------------------------------------------------------------------------------------------|
| **Input Type**        | Set to `fasta` when blasting defined sequences or `ab1` (default) for Sanger read input                                    |
| **Input Directory**   | The input directory should contain `.ab1` files, which are the raw Sanger sequencing files.      |
| **Filename Format**   | The script expects the filenames to follow a specific pattern: `<prefix1>_<prefix2>_<sampleID>_<orientation>_<optional_suffix>.ab1`  |
|             | `prefix1` / `prefix2`: Any string of characters (e.g., project name, date).                                             |
|             | `sampleID`: Unique identifier for the sample.                                                                |
|          | `orientation`: String containing 'F' for forward reads and 'R' for reverse reads.                                         |
|          | `optional_suffix`: Any string or strings following 3rd underscore                                         |

#### Example Filenames
- `user1_project1_sample1_27F_H02.ab1`
- `user1_project1_sample1_1492R_H02.ab1`


### Output

The script generates several output files based on the specified or default arguments:

| Output File             | Description                                                                            | Default Filename                      |
|-------------------------|----------------------------------------------------------------------------------------|---------------------------------------|
| BLAST Results           | Contains the results of the BLAST search for each sample.                              | `blast_results.csv`                   |
| Consensus FASTA         | FASTA file containing the consensus sequences generated from aligned reads.            | `consensus_sequences.fasta`           |
| Trimmed Consensus FASTA | FASTA file containing consensus sequences with large chunks of Ns trimmed.             | `trimmed_consensus_sequences.fasta`   |

The location of these output files can be controlled using the `--output_dir` argument. If `--output_dir` is not specified, the output files will be saved in the current directory.

#### BLAST Results Output Table

The following table describes the columns present in the BLAST results output CSV file:

| Column Name       | Description                                                                 |
|-------------------|-----------------------------------------------------------------------------|
| Query             | The identifier of the query sequence from the input FASTA file.             |
| Database          | The name of the BLAST database used for the search.                         |
| Query Length      | The length of the query sequence.                                           |
| Target Length     | The length of the target sequence.                                          |
| Alignment Length  | The length of the alignment between the query and target sequences.         |
| Score             | The score of the alignment. Higher scores indicate better alignments.       |
| E-value           | The expect value of the alignment. Lower values indicate more significant alignments. |
| Identity          | The number of identical matches in the alignment.                           |
| Percent Identity  | The percentage of identical matches in the alignment.                       |
| Query Coverage    | The percentage of the query sequence that is covered by the alignment.      |
| Mismatches        | The number of mismatched bases in the alignment.                            |
| Gap Openings      | The number of gaps in the alignment.                                        |
| Query Start       | The starting position of the alignment on the query sequence.               |
| Query End         | The ending position of the alignment on the query sequence.                 |
| Target Start      | The starting position of the alignment on the target sequence.              |
| Target End        | The ending position of the alignment on the target sequence.                |


### Databases

The script uses the following local BLAST databases for sequence alignment:

| Database Name      | Path                                                                                                       |  Description  |
|--------------------|------------------------------------------------------------------------------------------------------------|------|
| `16SMicrobial`     | `/usr2/people/kuehl/databases/16SMicrobial/16SMicrobial`                                                   | Sequences from NCBI database  |
| `CORAL_2024FEB`    | `/usr2/people/gtl/data/Isolate_Sanger_Checking/Collection_data/CORAL_2024FEB.fasta_db`                     |  Internal 16S sequences from Sanger data   |
| `Genome_16s`       | `/usr2/people/gtl/data/Isolate_Sanger_Checking/Collection_data/Genome.16s.fasta_db`                        |  Internal 16S seqeunces from genome assemblies   |


## Script Workflow

1. **Parse Arguments**: Parses the command line arguments.
2. **Get Files**: Reads .ab1 files from the specified input directory.
3. **Parse Filenames**: Identifies forward and reverse reads based on filenames.
4. **Process Sequences**:
    - Reads sequences from .ab1 files.
    - Aligns forward and reverse sequences to generate consensus sequences.
    - Optionally trims large chunks of Ns from consensus sequences.
5. **Save Sequences**: Writes consensus sequences to specified FASTA files.
6. **Perform BLAST**:
    - Executes BLAST searches against specified databases.
    - Parses BLAST results and saves to the output file.

## Functions

| Function | Description |
|----------|-------------|
| `parse_arguments()` | Parses command line arguments. |
| `get_sequence_from_ab1(file_path, verbose=False)` | Reads a sequence from an .ab1 file. |
| `align_sequences(seq1, seq2, verbose=False)` | Aligns two sequences and generates a consensus sequence. |
| `blast_sequence(sequence, blastdb=None, web_blast=False, verbose=False)` | Performs a BLAST search for a sequence. |
| `parse_blast_results(blast_record, verbose=False)` | Parses BLAST results. |
| `parse_filenames(file_list, verbose=False)` | Parses filenames to identify forward and reverse reads. |
| `get_files_from_directory(directory, suffix='.ab1', verbose=False)` | Gets .ab1 files from a directory. |
| `write_fasta(sequence, filename="query.fasta")` | Writes a sequence to a FASTA file. |


### Sequence Alignment

Alignment uses the `pairwise2.align.globalms` function from the Biopython library. It performs a global alignment between two sequences using specified match, mismatch, gap opening, and gap extension penalties.

#### Arguments

- `seq1`: The first sequence to align.
- `seq2`: The second sequence to align.
- `1`: The score for a match.
- `-1`: The score for a mismatch.
- `-2`: The penalty for opening a gap.
- `-0.5`: The penalty for extending a gap.

#### Description

- **Match Score (`1`)**: Adds 1 for each matching pair of bases.
- **Mismatch Penalty (`-1`)**: Subtracts 1 for each mismatched pair of bases.
- **Gap Opening Penalty (`-2`)**: Subtracts 2 when a gap is introduced in one of the sequences.
- **Gap Extension Penalty (`-0.5`)**: Subtracts 0.5 for each additional base in an existing gap.

These parameters are adjusted to be more lenient, allowing for a more flexible alignment that can tolerate mismatches and gaps better.

## Example 1

To process Sanger sequencing data in the directory `data` and save the results to `results.csv`:

```
python /usr2/people/anoonan/scripts/GitHub/sangerBLAST/sangerBLAST.py data --output results.csv --output_dir output_directory --consensus_fasta consensus_sequences.fasta --trimmed_fasta trimmed_sequences.fasta --verbose
```

This command will:

- Read `.ab1` files from the `data` directory.
- Generate consensus sequences.
- Perform local BLAST searches against the specified databases.
- Save the BLAST results to `output_directory/results.csv`.
- Save the consensus sequences to `output_directory/consensus_sequences.fasta` and `output_directory/trimmed_sequences.fasta`.

## Example 2

To process ASVs from FASTA file and save the results to `results.csv`:

```
python /usr2/people/anoonan/scripts/GitHub/sangerBLAST/sangerBLAST.py asv_seqs.fasta --output results.csv --output_dir output_directory --verbose
```

This command will:

- Import ASV sequences
- Perform local BLAST searches against the specified databases.
- Save the BLAST results to `output_directory/results.csv`.
- Save the consensus sequences to `output_directory/consensus_sequences.fasta` and `output_directory/trimmed_sequences.fasta`.


## License

This project is licensed under the MIT License.
