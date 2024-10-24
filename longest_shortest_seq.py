from Bio import SeqIO

# Path to your FASTA file
file_path = 'dna2.fasta'

# Dictionary to store sequences and their lengths
sequences = {}

# Read the file and get the lengths
for seq_record in SeqIO.parse(file_path, "fasta"):
    sequences[seq_record.id] = len(seq_record.seq)

# Find the longest and shortest sequence lengths
max_length = max(sequences.values())
min_length = min(sequences.values())

# Find all sequences with the max and min lengths
longest_seqs = [id for id, length in sequences.items() if length == max_length]
shortest_seqs = [id for id, length in sequences.items() if length == min_length]

print("Longest Sequence Length:", max_length)
print("Identifiers of Longest Sequences:", longest_seqs)
print("Shortest Sequence Length:", min_length)
print("Identifiers of Shortest Sequences:", shortest_seqs)
