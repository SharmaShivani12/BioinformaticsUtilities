count = 0
with open('dna2.fasta', 'r') as file:
    for line in file:
        if line.startswith('>'):
            count += 1
print(count)