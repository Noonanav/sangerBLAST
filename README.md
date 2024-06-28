# Sanger Sequencing Data Processor

This script processes Sanger sequencing data in a 2-step process: 

1. Paired forward and reverse reads are aligned using the `pairwise2.align.globalms` function to generate consensus sequences
2. A BLAST search against specified databases (either local and web-based) and top matches are reported in the output file.

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
sangerBLAST.py [-h] [--output OUTPUT] [--consensus_fasta CONSENSUS_FASTA] [--trimmed_fasta TRIMMED_FASTA] [--verbose] [--web_blast] [--trim_Ns] [--output_dir OUTPUT_DIR] input_directory

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
  --trim_Ns             Trim large chunks of Ns before BLAST
  --output_dir OUTPUT_DIR
                        Directory to place output files
```

### Arguments

| Argument            | Description                                                                    | Default                           |
|---------------------|--------------------------------------------------------------------------------|-----------------------------------|
| `--output`          | Output file for BLAST results                                                  | `blast_results.csv`               |
| `--output_dir`      | Directory to save output files                                                 | `.`                               |
| `--consensus_fasta` | Output file for consensus sequences in FASTA format                            | `consensus_sequences.fasta`       |
| `--trimmed_fasta`   | Output file for trimmed consensus sequences in FASTA format                    | `trimmed_consensus_sequences.fasta` |
| `--verbose`         | Enable verbose output                                                          |         `False`                          |
| `--web_blast`       | BLAST against NCBI database instead of local BLAST                             |     `False`                              |
| `--trim_Ns`         | Trim large chunks of Ns before BLAST                                           |      `False`                             |

### Input

The script assumes the following input format:

| Input Type            | Description                                                                                      |
|-----------------------|--------------------------------------------------------------------------------------------------|
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


### Databases

The script uses the following local BLAST databases for sequence alignment:

| Database Name      | Path                                                                                                       |
|--------------------|------------------------------------------------------------------------------------------------------------|
| `16SMicrobial`     | `/usr2/people/kuehl/databases/16SMicrobial/16SMicrobial`                                                   |
| `CORAL_2024FEB`    | `/usr2/people/gtl/data/Isolate_Sanger_Checking/Collection_data/CORAL_2024FEB.fasta_db`                     |
| `Genome_16s`       | `/usr2/people/gtl/data/Isolate_Sanger_Checking/Collection_data/Genome.16s.fasta_db`                        |


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

## Example

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


## License

This project is licensed under the MIT License.
