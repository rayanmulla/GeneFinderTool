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

def gene_finder(file_path):
    #Finding_genes_forward_or_reverse
    def find_genes_in_frame(sequence):
        def start_codon(i, sequence):
            return sequence[i:i+3] == 'ATG'
    
        def stop_codon(i, sequence):
            return sequence[i:i+3] in {'TAA', 'TAG', 'TGA'}
    
        #Storing_identified_gene_seqs
        gene_seqs = []
        #Screening_for_genes_with_start_and_stop_codon
        i = 0 
        while i < len(sequence) - 2:
            if start_codon(i, sequence):
                for j in range(i+3, len(sequence) -2, 3):
                    if stop_codon(j, sequence):
                        gene_seqs.append(sequence[i:j+3])
                        #Move_i_forward_to_skip_this_found_gene_and_look_next_start
                        i = j
                        break
            i += 1
        return gene_seqs
        
    #Reading_the_file_and_extracting_sequences
    def read_fasta(file_path):
        with open(file_path, 'r') as file:
            sequence = ""
            for line in file:
                if not line.startswith('>'):
                    sequence += line.strip()
        return sequence
    
    #Read_input
    dna_sequence = read_fasta(file_path)

    #Collect_potential_genes
    all_genes_seqs = []

    #Forward_analysis
    for frame in range(3):
        frame_sequence = dna_sequence[frame:]
        all_genes_seqs.extend(find_genes_in_frame(frame_sequence))

    #Reverse_analysis
    reverse_complement = str(Seq(dna_sequence).reverse_complement())
    for frame in range(3):
        frame_sequence = reverse_complement[frame:]
        all_genes_seqs.extend(find_genes_in_frame(frame_sequence))

    return all_genes_seqs
    
if __name__ == "__main__":
    input_file = "GCA_000006745.1_ASM674v1_genomic.fna"
    genes = gene_finder(input_file)
    if genes:
        print("Found gene seqs:")
        for gene in genes:
            print(gene)
    else:
        print("No genes found")