from Bio.Seq import Seq
from Bio import SeqIO
import sys

def translate_orfs(sequence, min_length=100):
    """Translate ORFs in all six reading frames to protein sequences,
    filtering ORFs based on minimum length provided in codons."""
    proteins = set()  # Use a set to store unique protein sequences

    def find_orfs(seq):
        """Find ORFs in a single sequence reading frame."""
        orfs = []  # To store ORFs as sequences
        start_codon = 'ATG'
        stop_codons = {'TAA', 'TAG', 'TGA'}

        i = 0
        while i < len(seq) - 2:
            if seq[i:i+3] == start_codon:
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
    # Get the file path from command line arguments
    file_path = sys.argv[1]

    # Default minimum ORF length is 100 codons, override if passed as argument
    min_length = int(sys.argv[2]) if len(sys.argv) > 2 else 100

    sequence = read_fasta(file_path)
    protein_strings = translate_orfs(sequence, min_length=min_length)
    
    # Print only the protein strings or indicate if none were found
    if protein_strings:
        for protein in protein_strings:
            print(protein)
    else:
        print("No protein sequences found.")
