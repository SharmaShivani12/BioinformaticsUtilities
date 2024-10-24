from Bio import SeqIO
from Bio.Seq import Seq

def find_orfs(seq, frame):
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    
    # Scan the sequence for ORFs in the specified frame
    for i in range(frame, len(seq), 3):
        codon = seq[i:i+3]
        if codon == start_codon:
            for j in range(i+3, len(seq), 3):
                if seq[j:j+3] in stop_codons:
                    orfs.append((i, seq[i:j+3]))
                    break
    return orfs

def analyze_sequences(file_path):
    longest_orf = None
    longest_orf_length = 0
    sequence_data = {}

    # Process each sequence in the FASTA file
    for record in SeqIO.parse(file_path, "fasta"):
        seq_id = record.id
        sequence_data[seq_id] = {
            'longest_orf_length': 0,
            'longest_orf_details': None
        }

        # Check all three forward frames
        for frame in range(3):
            orfs = find_orfs(str(record.seq), frame)
            for start, orf_seq in orfs:
                orf_length = len(orf_seq)
                if orf_length > sequence_data[seq_id]['longest_orf_length']:
                    sequence_data[seq_id]['longest_orf_length'] = orf_length
                    sequence_data[seq_id]['longest_orf_details'] = (start, orf_seq)
                
                if orf_length > longest_orf_length:
                    longest_orf_length = orf_length
                    longest_orf = (seq_id, start, orf_seq)

    return sequence_data, longest_orf

# Path to your FASTA file
file_path = 'dna2.fasta'
sequence_data, longest_orf = analyze_sequences(file_path)

print("Details of the longest ORFs in each sequence:")
for seq_id, details in sequence_data.items():
    print(f"Sequence ID: {seq_id}, Length: {details['longest_orf_length']}, Start: {details['longest_orf_details'][0] if details['longest_orf_details'] else 'N/A'}")

print("\nOverall longest ORF:")
print(f"Sequence ID: {longest_orf[0]}, Start: {longest_orf[1]}, ORF Length: {len(longest_orf[2])}")
