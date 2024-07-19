import os
import argparse
import pandas as pd
from collections import defaultdict
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
import time

blast_dbs = {
    '16SMicrobial': '/usr2/people/kuehl/databases/16SMicrobial/16SMicrobial',
    'CORAL_2024FEB': '/usr2/people/gtl/data/Isolate_Sanger_Checking/Collection_data/CORAL_2024FEB.fasta_db',
    'Genome_16s': '/usr2/people/gtl/data/Isolate_Sanger_Checking/Collection_data/Genome.16s.fasta_db'
}

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process Sanger sequencing data, align reads, generate consensus sequences, and perform BLAST.')
    parser.add_argument('input_directory', type=str, help='Directory containing .ab1 files')
    parser.add_argument('--output', type=str, default='blast_results.csv', help='Output file for BLAST results')
    parser.add_argument('--consensus_fasta', type=str, default='consensus_sequences.fasta', help='Output file for consensus sequences in FASTA format')
    parser.add_argument('--trimmed_fasta', type=str, default='trimmed_consensus_sequences.fasta', help='Output file for trimmed consensus sequences in FASTA format')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output', default=True)
    parser.add_argument('--web_blast', action='store_true', help='BLAST against NCBI database instead of local BLAST')
    parser.add_argument('--blast_Ns', action='store_true', help='BLAST concensus sequences with Ns')
    parser.add_argument('--output_dir', type=str, default='.', help='Directory to place output files')
    parser.add_argument('--input_type', choices=['ab1', 'fasta'], default='ab1', help='Type of input data (default: ab1)')
    return parser.parse_args()

def get_sequence_from_ab1(file_path, verbose=False):
    try:
        if verbose:
            print(f"\nReading sequence from {file_path}")
        record = SeqIO.read(file_path, 'abi')
        seq = record.seq
        
        # Optional: Filter low-quality bases if necessary
        if "phred_quality" in record.letter_annotations:
            if verbose:
                print("Filtering low-quality bases...")
            quality_scores = record.letter_annotations["phred_quality"]
            filtered_seq = ''.join(
                base if quality >= 20 else 'N' for base, quality in zip(seq, quality_scores)
            )
            seq = Seq(filtered_seq)
        
        return seq
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def align_sequences(seq1, seq2, verbose=False):
    # Reverse complement the reverse sequence
    seq2_rc = seq2.reverse_complement()
    if verbose:
        print("\nAligning sequences...")
    # Adjusting alignment parameters to be more lenient
    alignments = pairwise2.align.globalms(seq1, seq2_rc, 1, -1, -2, -0.5)  # Match: 1, Mismatch: -1, Gap open: -2, Gap extend: -0.5
    if alignments:
        top_alignment = alignments[0]
        aligned_seq1, aligned_seq2 = top_alignment.seqA, top_alignment.seqB
        consensus_seq = ''.join(
            base1 if base1 == base2 else 'N' if base1 != 'N' and base2 != 'N' else base1 if base1 != 'N' else base2
            for base1, base2 in zip(aligned_seq1, aligned_seq2)
        )
        if verbose:
            print(f"Consensus sequence: {consensus_seq}")
        return Seq(consensus_seq)  # Return as a Seq object
    else:
        print("No alignments found.")
        return None

def trim_large_chunks_of_ns(sequence, min_ns_stretch=10):
    import re

    # Find all stretches of min_ns_stretch or more Ns
    stretches = [(m.start(), m.end()) for m in re.finditer(r'N{' + str(min_ns_stretch) + r',}', sequence)]
    
    if not stretches:
        return sequence  # No large chunks of Ns found

    # Calculate gaps between the end of one stretch and the start of the next
    gaps = [(stretches[i][1], stretches[i+1][0]) for i in range(len(stretches) - 1)]
    
    # Identify the largest gap
    largest_gap = max(gaps, key=lambda x: x[1] - x[0])
    
    # Extract the trimmed sequence around the largest gap
    trimmed_sequence = sequence[largest_gap[0]:largest_gap[1]]
    return trimmed_sequence

def blast_sequence(sequence, output_dir, blastdb=None, web_blast=False, verbose=False):
    if web_blast:
        try:
            if verbose:
                print("\nPerforming BLAST search against NCBI...")
            result_handle = NCBIWWW.qblast("blastn", "refseq_rna", sequence)
            if verbose:
                print("BLAST search completed, parsing results...")
            blast_record = NCBIXML.read(result_handle)
            return blast_record
        except Exception as e:
            print(f"Error during BLAST: {e}")
            return None
    else:
        if verbose:
            print("\nPerforming local BLAST search...")
        query_file_path = os.path.join(output_dir, "query.fasta")
        output_file = os.path.join(output_dir, "blast_output.xml")
        try:
            write_fasta(sequence, query_file_path)  # Use the new function to write the FASTA file
            blastn_cline = NcbiblastnCommandline(query=query_file_path, db=blastdb, evalue=0.001, outfmt=5, out=output_file)
            stdout, stderr = blastn_cline()
            blast_record = None
            if os.path.exists(output_file):
                with open(output_file) as result_handle:
                    blast_record = NCBIXML.read(result_handle)
                os.remove(output_file)
            os.remove(query_file_path)
            return blast_record
        except ApplicationError as e:
            print(f"ApplicationError during BLAST: {e}")
            if os.path.exists(query_file_path):
                os.remove(query_file_path)
            if os.path.exists(output_file):
                os.remove(output_file)
            return None
        except Exception as e:
            print(f"Error during BLAST: {e}")
            if os.path.exists(query_file_path):
                os.remove(query_file_path)
            if os.path.exists(output_file):
                os.remove(output_file)
            return None

def parse_blast_results(blast_record, verbose=False):
    top_hits = []
    if blast_record:
        if verbose:
            print("\nParsing BLAST results...")
        hits_count = 0
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hits_count >= 5:
                    break
                # Calculate percent identity and query coverage
                percent_identity = (hsp.identities / hsp.align_length) * 100
                query_coverage = ((hsp.query_end - hsp.query_start + 1) / blast_record.query_length) * 100

                # Calculate alignment length manually
                alignment_length = hsp.query_end - hsp.query_start + 1

                hit = {
                    'Target': alignment.title,
                    'Query Length': blast_record.query_length,
                    'Target Length': alignment.length,
                    'Alignment Length': alignment_length,
                    'Score': hsp.score,
                    'E-value': hsp.expect,
                    'Identity': hsp.identities,
                    'Percent Identity': percent_identity,
                    'Query Coverage': query_coverage,
                    'Mismatches': hsp.align_length - hsp.identities,
                    'Gap Openings': hsp.gaps,
                    'Query Start': hsp.query_start,
                    'Query End': hsp.query_end,
                    'Target Start': hsp.sbjct_start,
                    'Target End': hsp.sbjct_end
                }
                top_hits.append(hit)
                hits_count += 1
                if verbose:
                    print(f"\nTarget: {hit['Target']}")
                    print(f"Query Length: {hit['Query Length']}")
                    print(f"Target Length: {hit['Target Length']}")
                    print(f"Alignment Length: {hit['Alignment Length']}")
                    print(f"Score: {hit['Score']}")
                    print(f"E-value: {hit['E-value']}")
                    print(f"Identity: {hit['Identity']}")
                    print(f"Percent Identity: {hit['Percent Identity']:.2f}%")
                    print(f"Query Coverage: {hit['Query Coverage']:.2f}%")
                    print(f"Mismatches: {hit['Mismatches']}")
                    print(f"Gap Openings: {hit['Gap Openings']}")
                    print(f"Query Start: {hit['Query Start']}")
                    print(f"Query End: {hit['Query End']}")
                    print(f"Target Start: {hit['Target Start']}")
                    print(f"Target End: {hit['Target End']}")
                    print('-' * 60)
            if hits_count >= 5:
                break
    return top_hits

def parse_filenames(file_list, verbose=False):
    sequences = defaultdict(lambda: {'fwd': None, 'rev': None})
    
    for file_name in file_list:
        parts = file_name.split('_')
        if len(parts) >= 4:
            sample_id = parts[2]
            orientation = parts[3]
            if 'F' in orientation:
                sequences[sample_id]['fwd'] = file_name
                if verbose:
                    print(f"Forward read identified: {file_name}")
            elif 'R' in orientation:
                sequences[sample_id]['rev'] = file_name
                if verbose:
                    print(f"Reverse read identified: {file_name}")
        else:
            print(f"Invalid filename format: {file_name}")
    
    return sequences

def get_files_from_directory(directory, suffix='.ab1', verbose=False):
    files = [file for file in os.listdir(directory) if file.endswith(suffix)]
    if verbose:
        print(f"Found {len(files)} .ab1 files in {directory}")
    return files

def write_fasta(sequence, filename="query.fasta"):
    with open(filename, "w") as f:
        f.write(f">query\n{sequence}\n")

def main():
    print("\n########################################")
    print("### Sanger Sequencing Data Processor ###")
    print("########################################\n")
    
    args = parse_arguments()
    input_directory = args.input_directory
    output_dir = args.output_dir
    output_file = os.path.join(args.output_dir, args.output)
    consensus_fasta = os.path.join(args.output_dir, args.consensus_fasta)
    trimmed_fasta = os.path.join(args.output_dir, args.trimmed_fasta)
    web_blast = args.web_blast
    blast_Ns = args.blast_Ns
    verbose = args.verbose
    input_type = args.input_type

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if verbose:
        print("\n========== Configuration ==========")
        print(f"Input directory: {input_directory}")
        print(f"Output file: {output_file}")
        print(f"Consensus FASTA file: {consensus_fasta}")
        print(f"Trimmed FASTA file: {trimmed_fasta}")
        print(f"BLAST against NCBI: {web_blast}")
        print(f"BLAST consensus sequence with Ns: {blast_Ns}")
        print(f"Input type: {input_type}")
        print("====================================\n")

    if input_type == 'ab1':
        file_list = [file for file in os.listdir(input_directory) if file.endswith('.ab1')]
        sequences = parse_filenames(file_list, verbose=verbose)

        print("\nProcessing", len(sequences), "samples...")

        consensus_sequences = {}
        consensus_records = []
        trimmed_records = []
        for sample_id, seqs in sequences.items():
            fwd_file = seqs['fwd']
            rev_file = seqs['rev']
            if fwd_file and rev_file:
                fwd_seq = get_sequence_from_ab1(os.path.join(input_directory, fwd_file), verbose=verbose)
                rev_seq = get_sequence_from_ab1(os.path.join(input_directory, rev_file), verbose=verbose)
                if fwd_seq and rev_seq:
                    consensus_seq = align_sequences(fwd_seq, rev_seq, verbose=verbose)
                    if consensus_seq:
                        consensus_length = len(consensus_seq)
                        consensus_records.append(SeqRecord(consensus_seq, id=sample_id, description=f"Consensus sequence length={consensus_length}"))
                        trimmed_seq = trim_large_chunks_of_ns(str(consensus_seq))
                        trimmed_length = len(trimmed_seq)
                        trimmed_records.append(SeqRecord(Seq(trimmed_seq), id=sample_id, description=f"Trimmed consensus sequence length={trimmed_length}"))
                        if blast_Ns:
                            consensus_sequences[sample_id] = consensus_seq
                        else:
                            consensus_sequences[sample_id] = Seq(trimmed_seq)

        # Write consensus sequences to FASTA files
        if consensus_records:
            with open(consensus_fasta, 'w') as fasta_output:
                SeqIO.write(consensus_records, fasta_output, 'fasta')
            if verbose:
                print(f"\nConsensus sequences saved to {consensus_fasta}")

        if trimmed_records:
            with open(trimmed_fasta, 'w') as fasta_output:
                SeqIO.write(trimmed_records, fasta_output, 'fasta')
            if verbose:
                print(f"\nTrimmed consensus sequences saved to {trimmed_fasta}")

    elif input_type == 'fasta':
        # Assume input_directory is a fasta file containing multiple sequences
        sequences = SeqIO.parse(input_directory, 'fasta')
        consensus_sequences = {seq.id: seq.seq for seq in sequences}
        # No consensus sequence or trimming step needed

    full_blast_results = pd.DataFrame()
    for db_name, blastdb in blast_dbs.items():
        if verbose:
            print(f"\nPerforming BLAST search against {db_name}")
        blast_results = []
        for sample_id, sequence in consensus_sequences.items():
            if verbose:
                print(f"\nPerforming BLAST for sample: {sample_id}")
            blast_record = blast_sequence(sequence, output_dir, blastdb=blastdb, web_blast=web_blast, verbose=verbose)
            if blast_record:
                top_hits = parse_blast_results(blast_record, verbose=verbose)
                for hit in top_hits:
                    hit['Query'] = sample_id
                    hit['Database'] = db_name
                blast_results.extend(top_hits)
                time.sleep(5)  # Add delay to prevent being blocked by NCBI servers
            else:
                if verbose:
                    print(f"\nFailed to get BLAST results for sample: {sample_id}")

        df = pd.DataFrame(blast_results)

        cols = ['Query'] + [col for col in df.columns if col != 'Query']
        df = df[cols]
        full_blast_results = pd.concat([full_blast_results, df])

    full_blast_results = full_blast_results.sort_values(by=['Query', 'Database', 'Score'], ascending=[True, True, False])

    full_blast_results.to_csv(output_file, index=False)
    if verbose:
        print(f"\nResults saved to {output_file}")

    if verbose:
        print("\n========== Summary ==========")
        print(f"Total sequences processed: {len(consensus_sequences)}")
        print(f"Total BLAST hits found: {len(full_blast_results)}")
        print("=============================\n")

if __name__ == '__main__':
    main()
