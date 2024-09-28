#Gene Finder Tool
#/home/mullara/genomes/class_genomes
#/home/mullara/GeneFinder/class_genomes/GCA_000006745.1/GCA_000006745.1_ASM674v1_genomic.fna
#LLM used ChatGPT-4o
#Prompt: make me a python tool for screening fasta files for genes that start with a start codon ('ATG') 
#and ends with one of the stop codons ('TAA', 'TAG', and 'TGA'). the tool inputs a .fna file 
#and outputs any region between the start codon and one of the end codons


#gene_finder()
#gene_finder_including_reverse()

from Bio.Seq import Seq
from Bio import SeqIO

def translate_orfs(sequence):
    """Translate ORFs in all six reading frames to protein sequences."""
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

    # Analyze forward strands
    for frame in range(3):
        frame_seq = sequence[frame:]
        orfs = find_orfs(frame_seq)
        for orf in orfs:
            orf_seq = Seq(orf)
            protein = orf_seq.translate(to_stop=True)
            proteins.add(str(protein))

    # Analyze reverse complement strands
    reverse_complement = str(Seq(sequence).reverse_complement())
    for frame in range(3):
        frame_seq = reverse_complement[frame:]
        orfs = find_orfs(frame_seq)
        for orf in orfs:
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
    input_file = "rosalind_orf_2.fna"  # Replace with your FASTA file name
    sequence = read_fasta(input_file)
    protein_strings = translate_orfs(sequence)
    
    # Print only the protein strings or indicate if none were found
    if protein_strings:
        for protein in protein_strings:
            print(protein)
    else:
        print("No protein sequences found.")


