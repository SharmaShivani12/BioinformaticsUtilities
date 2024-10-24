from Bio import SeqIO
from collections import defaultdict

def find_repeats(file_path, repeat_length=12):
    # Dictionary to store all repeats and their frequencies
    repeat_counts = defaultdict(int)

    # Read sequences from the FASTA file
    for record in SeqIO.parse(file_path, "fasta"):
        sequence = str(record.seq)
        # Use a sliding window to find all repeats of the given length
        for i in range(len(sequence) - repeat_length + 1):
            repeat = sequence[i:i + repeat_length]
            repeat_counts[repeat] += 1

    # Find the maximum frequency of any repeat
    max_frequency = max(repeat_counts.values())
    # Find all repeats that occur 'max_frequency' times
    max_repeats = [repeat for repeat, count in repeat_counts.items() if count == max_frequency]

    return max_frequency, max_repeats

# Path to your FASTA file
file_path = 'dna2.fasta'

# Find repeats of length 12 and determine the most frequent ones
max_frequency, max_repeats = find_repeats(file_path, 6)

print("Maximum frequency of any 12-base repeat:", max_frequency)
print("Number of different 12-base sequences that occur this maximum number of times:", len(max_repeats))
