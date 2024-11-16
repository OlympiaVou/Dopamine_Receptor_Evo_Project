'''
This script reduces the fasta sequences of a multi-fasta file.
INPUT:
1. A .txt file with the names of the species that the user wants
to include in the output file
2. A multifasta file
3. The name of the output file

OUTPUT: A fasta file only with the organisms of the .txt file
'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species", type=str)        # Save the .txt file
parser.add_argument("fasta", type=str)          # Save the input fasta file
parser.add_argument("output_file", type=str)    # Save the output fasta file
args = parser.parse_args()

input = args.fasta				# Save the input file name as input
out = args.output_file 				# Save the output file name as out
species_names = []                              # Initialize the list of species
with open(args.species,'r') as file:            # Open and read the .txt file
    for line in file:                           # For each line in the file
        species_names.append(line[:-1])         # Append the line in the list of species


from Bio import SeqIO # Library to process fasta files

records = list(SeqIO.parse(input, 'fasta')) # Read and save all the fasta entries of the input file

with open(out, 'w') as output:              # Open the output file
    for sp in species_names:                # For each species in the list of species
        for record in records:              # Check all the fasta entries of the input file
            if sp in record.description:    # If the species exists in the species list
                output.write(f'>{record.description}\n{record.seq}\n')  # Write the entry in the output file

