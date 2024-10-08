from Bio.Seq import Seq
from Bio import SeqIO
import sys
import os

def translate_orfs(sequence, min_length=100, rbs_distance=20):
    """Translate ORFs in all six reading frames to protein sequences,
    filtering ORFs based on minimum length and ribosome binding site presence."""
    proteins = set()  # Use a set to store unique protein sequences

    def find_orfs(seq):
        """Find ORFs in a single sequence reading frame."""
        orfs = []  # To store ORFs as sequences
        start_codon = 'ATG'
        stop_codons = {'TAA', 'TAG', 'TGA'}
        rbs_sequence = 'AGGAGG'  # Example Shine-Dalgarno sequence

        i = 0
        while i < len(seq) - 2:
            if seq[i:i+3] == start_codon:
                # Check for RBS presence upstream of start codon
                start_rbs_check = max(0, i - rbs_distance)
                upstream_sequence = seq[start_rbs_check:i]
                
                if rbs_sequence not in upstream_sequence:
                    i += 1
                    continue

                # Search for stop codon
                for j in range(i+3, len(seq) - 2, 3):
                    codon = seq[j:j+3]
                    if codon in stop_codons:
                        orf_seq = seq[i:j+3]
                        orfs.append(orf_seq)
                        break
            i += 1
        return orfs

    def valid_orf(orf_seq):
        """Check if the ORF sequence has the minimum number of codons."""
        return len(orf_seq) / 3 >= min_length

    # Analyze forward strands
    for frame in range(3):
        frame_seq = sequence[frame:]
        orfs = find_orfs(frame_seq)
        for orf in orfs:
            if valid_orf(orf):
                orf_seq = Seq(orf)
                protein = orf_seq.translate(to_stop=True)
                proteins.add(str(protein))

    # Analyze reverse complement strands
    reverse_complement = str(Seq(sequence).reverse_complement())
    for frame in range(3):
        frame_seq = reverse_complement[frame:]
        orfs = find_orfs(frame_seq)
        for orf in orfs:
            if valid_orf(orf):
                orf_seq = Seq(orf)
                protein = orf_seq.translate(to_stop=True)
                proteins.add(str(protein))

    return proteins

def read_fasta(file_path):
    """Read sequence from FASTA file."""
    with open(file_path, 'r') as file:
        for record in SeqIO.parse(file, "fasta"):
            return str(record.seq)

if __name__ == "__main__":
    # Get the file path and optional parameters from command line arguments
    file_path = sys.argv[1]
    output_directory = sys.argv[2]
    min_length = int(sys.argv[3]) if len(sys.argv) > 3 else 100
    rbs_distance = int(sys.argv[4]) if len(sys.argv) > 4 else 20

    sequence = read_fasta(file_path)
    protein_strings = translate_orfs(sequence, min_length=min_length, rbs_distance=rbs_distance)
    
    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)
    
    # Determine the output file path
    base_name = os.path.splitext(os.path.basename(file_path))[0]  # Remove the directory and file extension
    output_file_path = os.path.join(output_directory, base_name + "_AA_seq.txt")  # Create the full output file path

    # Write protein strings to the output file
    with open(output_file_path, 'w') as output_file:
        if protein_strings:
            for protein in protein_strings:
                output_file.write(protein + '\n')
        else:
            output_file.write("No protein sequences found.")

